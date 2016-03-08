/**
 * @file src/Zuker.cpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2016 Youri Hoogstrate
 *
 * This file is part of segmentation-fold and originally taken from
 * yh-kt-fold.
 *
 * segmentation-fold is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * segmentation-fold is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 * </PRE>
 */


#include "Pair.hpp"
#include "Region.hpp"
#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "Position.hpp"
#include "PairingPlus.hpp"
#include "SubSequence.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"
#include "SegmentLoop.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"
#include "SegmentLoopTree.hpp"

#include "ScoringMatrix.hpp"
#include "Settings.hpp"
#include "DotBracket.hpp"
#include "ReadData.hpp"

#include "Zuker.hpp"



/**
 * @brief Constructs /initializes the Zuker class: include parameters and init an empty dotbracket output.
 *
 * @todo move this to this->init(); and run this->init(); or rename it to this->reset();
 */
Zuker::Zuker(Settings &arg_settings, Sequence &arg_sequence, ReadData &arg_thermodynamics) :
	GibbsFreeEnergy(arg_sequence, arg_thermodynamics),
	settings(arg_settings),

	vij(arg_sequence.size(), N_INFINITY),
	wij(arg_sequence.size(), 0.0),
	wmij(arg_sequence.size(), N_INFINITY),

	tij_v(arg_sequence.size(), {Pair(UNBOUND, UNBOUND), V_MATRIX}),				//@todo use N instead of 0? >> if so, set UNBOUND to N  + 1 or so
	  tij_w(arg_sequence.size(), {Pair(UNBOUND, UNBOUND), W_MATRIX}),				//@todo use N instead of 0? >> if so, set UNBOUND to N  + 1 or so
	  tij_wm(arg_sequence.size(), {Pair(UNBOUND, UNBOUND), WM_MATRIX}),			//@todo use N instead of 0? >> if so, set UNBOUND to N  + 1 or so

	  sij(arg_sequence.size(), nullptr)
{
#if DEBUG
	this->tij_v.fill({{NOT_YET_CALCULATED, NOT_YET_CALCULATED}, V_MATRIX});
	this->tij_w.fill({{NOT_YET_CALCULATED, NOT_YET_CALCULATED}, W_MATRIX});
#endif //DEBUG
	
	this->sequence_begin = this->sequence.data.begin();
	this->traceback_stacktop = -1;
}



/**
 * @brief Calculate the Gibbs free energy; fill the matrices Vij and Wij.
 *
 * @section DESCRIPTION a new format, having been standardized in 2008 based
 * The Vij, Wij and Mij matrix are filled in this function. This is done
 * by walking over the diagonals of the matrix to avoid a recursive
 * explosion. Each position in a diagnoal can be calculated independently
 * from each other.
 *
 * @todo Return: energy at i,j
 */
float Zuker::energy(void)
{
	unsigned int i = 0;
	unsigned int k = 0;
	
	for(k = 1; k < this->sequence.size(); k++)
	{
		// Paralelization / threading: "still reachable" memory error seems normal (http://stackoverflow.com/questions/6973489/valgrind-and-openmp-still-reachable-and-possibly-lost-is-that-bad).. -num_threads can be defined here
		// http://people.cs.pitt.edu/~melhem/courses/xx45p/OpenMp.pdf
		// For each element in a diagonal, it can be calculated independent of the others - if num_threads = 0, it will take all possible threads
		
		#pragma omp parallel for num_threads(this->settings.num_threads) private(i)
		for(i = 0; i < this->sequence.size() - k; i++)
		{
			Pair pair(i, i + k);// j = i + k
			this->w(pair);//this->w(i, i + k);
		}
	}
	
	Pair pair = Pair(0, this->sequence.size() - 1);
	return this->wij.get(pair);
}



/**
 * @brief Vij Function - energy if (i,j) pair, otherwise return infinity
 *
 * @param p1 A pair of positions refering to Nucleotide positions in the sequence, where pi.first < p1.second
 *
 * @return amount of Gibbs free energy provided for folding nucleotide i with j assuming i and j are paired
 */
float Zuker::v(Pair &p1, PairingPlus &p1p)
{
#if DEBUG
	if(this->tij_v.get(p1).target.first != (unsigned int) NOT_YET_CALCULATED)
	{
		throw std::invalid_argument("Zuker::v(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + "): redundant calculation, please request values from the ScoringMatrix directly");
	}
	
	if(!p1p.is_canonical() || (p1p.size) < this->settings.minimal_hairpin_length)
	{
		throw std::invalid_argument("Zuker::v(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + "): this pair should never be checked within this function because it's energy is infinity by definition");
	}
#endif //DEBUG
	
	float energy, tmp, tmp_k;
	
	Segment *tmp_segment;
	SegmentLoop *tmp_segmentloop;
	SegmentTraceback  *tmp_segmenttraceback = nullptr;
	
	energy = this->get_hairpin_loop_element(p1);						// Hairpin element
	traceback_jump tmp_tij = traceback_jump {{UNBOUND, UNBOUND}, V_MATRIX};
	
	
	// Stacking element
	Pair p2 = Pair(p1.first + 1, p1.second - 1);
	float v_ij_jp = this->vij.get(p2);
	tmp = this->get_stacking_pair_element(p1) + v_ij_jp;
	if(tmp < energy)
	{
		energy = tmp;
		tmp_tij.target = p2;
	}
	
	// SegmentLoop element
	SubSequence ps1 = this->sequence.ssubseq(p2.first, p2.second); ///@todo use Pair()
	tmp_segmentloop = this->thermodynamics.segmentloops.search(ps1);
	if(tmp_segmentloop != nullptr)
	{
		tmp_k = tmp_segmentloop->gibbs_free_energy + this->get_stacking_pair_without_surrounding(p1p);
		if(tmp_k < energy)
		{
			energy = tmp_k;
			tmp_segmenttraceback = &tmp_segmentloop->traceback;
		}
	}
	
	// L-bulge
	p2.second = p1.second - 1;
	for(p2.first = p1.first + 2; p2.first < p2.second; p2.first++)
	{
		PairingPlus pairing2 = PairingPlus(this->sequence_begin + p2.first, this->sequence_begin + p2.second);
		if(pairing2.is_canonical())
		{
			v_ij_jp = this->vij.get(p2);
			Region region = Region {p1, p2};
			
			tmp = this->get_bulge_loop_element(region) + v_ij_jp;
			if(tmp < energy)
			{
				energy = tmp;
				tmp_tij.target = p2;
				//tmp_segmenttraceback = nullptr;
			}
		}
	}
	
	// R-bulge
	p2.first = p1.first + 1;
	for(p2.second = p1.second - 2; p2.second > p2.first; p2.second--)
	{
		PairingPlus pairing2 = PairingPlus(this->sequence_begin + p2.first, this->sequence_begin + p2.second);
		if(pairing2.is_canonical())
		{
			v_ij_jp = this->vij.get(p2);
			Region region = Region {p1, p2};
			
			tmp = this->get_bulge_loop_element(region) + v_ij_jp;
			if(tmp < energy)
			{
				energy = tmp;
				tmp_tij.target = p2;
				//tmp_segmenttraceback = nullptr;
			}
		}
	}
	
	// Interior or Segment
	for(p2.first = p1.first + 2; p2.first < p1.second; p2.first++)
	{
		for(p2.second = p1.second - 2; p2.second > p2.first; p2.second--)
		{
			PairingPlus pairing2 = PairingPlus(this->sequence_begin + p2.first, this->sequence_begin + p2.second);
			if(pairing2.is_canonical())									// The following structure elements must be enclosed by pairings on both sides
			{
				v_ij_jp = this->vij.get(p2);
				Region region = Region {p1, p2};
				
				tmp = this->get_interior_loop_element(region) + v_ij_jp;
				if(tmp < energy)
				{
					energy = tmp;
					
					tmp_tij.target = p2;
					tmp_segmenttraceback = nullptr;
					//tmp_tij.target_matrix = V_MATRIX;
				}
				
				// Find segments:
				SubSequence pp1 = this->sequence.ssubseq(p1.first + 1, p2.first - 1);
				SubSequence pp2 = this->sequence.ssubseq(p2.second + 1, p1.second - 1);
				tmp_segment = this->thermodynamics.segments.search(pp1 , pp2);
				
				if(tmp_segment != nullptr)
				{
					tmp_k = tmp_segment->gibbs_free_energy + this->get_stacking_pair_without_surrounding(p1p) + v_ij_jp;
					
					if(tmp_k < energy)
					{
						energy = tmp_k;
						
						tmp_tij.target = p2;
						tmp_segmenttraceback = &tmp_segment->traceback;
						//tmp_tij.target_matrix = V_MATRIX;
					}
				}
			}
		}
	}
	
	
	//this->settings.minimal_hairpin_length)
	/*   CAA
	 * CC \ A
	 * ||   G   <- k (p2.first)
	 * CC  /G
	 *   C  A
	 *    AA
	 *
	 */
	// Multi-loop - comparison part of the loop looks weird, but otherwise it may become negative while it's unsigned
	for(
		p2.first = p1.first + 2 + this->settings.minimal_hairpin_length;
		p2.first + 2 + this->settings.minimal_hairpin_length < p1.second;
		p2.first++)
	{
#if DEBUG
		///@todo use minimal hairpin size as well - be aware of unsigned int
		if(p1.first + 1 >= p2.first || p2.first >= p1.second - 2)
		{
			throw std::invalid_argument("Zuker::v(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + ") -> m(" + std::to_string(p1.first + 1) + ", " + std::to_string(p2.first) + ", " + std::to_string(p2.first + 1) + ", " + std::to_string(p1.second - 1) + "): incorrect loop size");
		}
#endif //DEBUG
		Pair p3 = Pair(p1.first + 1, p2.first);
		Pair p4 = Pair(p2.first + 1, p1.second - 1);
		tmp = this->wmij.get(p3) + this->wmij.get(p4);
		
		if(tmp < energy)
		{
			energy = tmp;
			
			tmp_tij.target = {p2.first, p2.first};
			tmp_tij.target_matrix = WM_MATRIX;
			//tmp_segmenttraceback = nullptr;
		}
	}
	
#if DEBUG
	if(tmp_segmenttraceback != nullptr && tmp_tij.target_matrix == WM_MATRIX)
	{
		throw std::invalid_argument("Zuker::v(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + ") impossible: bifucation && segment");
	}
#endif //DEBUG
	
	this->tij_v.set(p1, tmp_tij);
	if(tmp_segmenttraceback != nullptr && tmp_tij.target_matrix != WM_MATRIX)
	{
		this->sij.set(p1, tmp_segmenttraceback);
	}
	
	this->vij.set(p1, energy);
	return energy;
}



/**
 * @brief Wij Function - provided Gibbs free energy for sequence i...j
 *
 * @param p1 A pair of positions refering to Nucleotide positions in the sequence, where pi.first < p1.second
 *
 * @return amount of Gibbs free energy provided for folding nucleotide i with j
 */
float Zuker::w(Pair &p1)
{
#if DEBUG
	if(p1.first >= p1.second)
	{
		throw std::invalid_argument("Zuker::w(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + "): out of bound");
	}
	
	if(this->tij_w.get(p1).target.first != (unsigned int) NOT_YET_CALCULATED)
	{
		throw std::invalid_argument("Zuker::w(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + "): redundant calculation, please request values from the ScoringMatrix directly");
	}
#endif //DEBUG
	
	float energy;
	signed int tmp_pij, tmp_qij;
	unsigned int n = (p1.second - p1.first);
	
	
	energy = 0.0;
	tmp_pij = UNBOUND;
	PairingPlus p1p = PairingPlus(this->sequence_begin + p1.first, this->sequence_begin + p1.second);
	
	if(n <= this->settings.minimal_hairpin_length)
	{
		this->vij.set(p1, N_INFINITY);
	}
	else
	{
		float tmp;
		
		if(!p1p.is_canonical())
		{
			this->vij.set(p1, N_INFINITY);
		}
		else
		{
			tmp = this->v(p1, p1p);
			if(tmp < energy)
			{
				energy = tmp;
				tmp_pij = BOUND;
			}
		}
		
		
		/*
		The following bifurcations have to be tested:
		
		CC (n = 1) -> none
		||
		
		CCC
		 (|  (pre-iter 1)
		|)   (pre-iter 2)
		
		CCCAGGA
		 (    |  (pre-iter 1)
		|    )   (pre-iter 2)
		|)(   |  (iter 1; k=+1)
		| )(  |  (iter 2; k=+2)
		|  )( |  (iter 3; k=+3)
		|   )(|  (iter 4; k=+4)
		*/
		if(n >= 2)// && tmp_pij != BOUND)// if it is bound, use Vij
		{
			Pair p2, p3;
			
			// pre-iter 1
			p3 = Pair(p1.first + 1, p1.second);
			tmp = this->wij.get(p3);
			if(tmp < energy)
			{
				energy = tmp;
				
				tmp_pij = (signed int) p3.first;
				tmp_qij = (signed int) p3.second;
			}
			
			// pre-iter 2
			p2 = Pair(p1.first, p1.second - 1);
			tmp = this->wij.get(p2);
			if(tmp < energy)
			{
				energy = tmp;
				tmp_pij = (signed int) p2.first;
				tmp_qij = (signed int) p2.second;
			}
			
			// remaining iterations
			for(unsigned int k = p1.first + 1; k < p1.second - 1; k++)
			{
				p2 = Pair(p1.first, k);
				p3 = Pair(k + 1, p1.second);
				
				tmp = this->wij.get(p2) + this->wij.get(p3);
				
				if(tmp < energy)
				{
					// Can also be done by checking whether pathmatrix_corrected_from > 0 ? >> and only store those positions in a tree instead of an entire matrix
					energy = tmp;
					tmp_pij = (signed int) k;///@todo make tmp_pij unsigned by setting BOUND and UNBOUND to MAX_VAL(size_t)-1 and MAX_VAL(size_t)-2
					tmp_qij = (signed int) k;///@todo make tmp_pij unsigned by setting BOUND and UNBOUND to MAX_VAL(size_t)-1 and MAX_VAL(size_t)-2
				}
			}
		}
	}
	
	// Calculate it, must be after v()
	this->wm(p1, p1p);
	this->wij.set(p1, energy);
	
	if(tmp_pij == BOUND)// no oneliner possible
	{
		this->tij_w.set(p1, {p1, V_MATRIX});
	}
	else
	{
		this->tij_w.set(p1, {{tmp_pij, tmp_qij}, W_MATRIX});
	}
	
	return energy;
}



/**
 * @brief Decomposes a multiloop
 */
float Zuker::wm(Pair &p1, PairingPlus &p1p)
{
#if DEBUG
	if(p1.first >= p1.second)
	{
		throw std::invalid_argument("Zuker::wm(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + "): out of bound");
	}
#endif //DEBUG
	
	unsigned int k;
	float tmp;
	float energy;
	energy = this->vij.get(p1);
	
	traceback_jump tmp_tij = {p1, V_MATRIX};
	Pair p2, p3;
	
	for(k = p1.first + 1; k < p1.second; k++)
	{
		p2 = Pair(p1.first, k);
		p3 = Pair(k + 1, p1.second);
		tmp = this->wmij.get(p2) + this->wmij.get(p3);
		
		if(tmp < energy)
		{
			// TB to p2 and p3, as bifurcation
			energy = tmp;
			tmp_tij.target = {k, k};
			tmp_tij.target_matrix = WM_MATRIX; // go back to the v matrix - for both
		}
	}
	
	p2 = Pair(p1.first + 1, p1.second);
	tmp = this->wmij.get(p2);
	if(tmp < energy)
	{
		// TB to p1l, not as bifurcation
		energy = tmp;
		tmp_tij.target = p2;
		tmp_tij.target_matrix = WM_MATRIX;//stay in wm
	}
	
	p2 = Pair(p1.first , p1.second - 1);
	tmp = this->wmij.get(p2);
	if(tmp < energy)
	{
		// TB to p1r, not as bifurcation
		energy = tmp;
		tmp_tij.target = p2;
		tmp_tij.target_matrix = WM_MATRIX;//stay in wm
	}
	
	this->tij_wm.set(p1, tmp_tij);
	this->wmij.set(p1, energy);
	
	return energy;
}



/**
 * @brief The traceback algorithm, finds the optimal path through the matrices.
 *
 * @section DESCRIPTION
 * This function finds the path back. It will choose between
 * the route provided by V or W scoring.
 *
 * The usage of an additional push and pop system is essential because
 * traces can split up because of forks. Otherwise recursion was
 * essential.
 */
void Zuker::traceback(void)
{
	this->folded_segments = 0;
	
	SegmentTraceback *independent_segment_traceback;
	Pair pair_tmp;
	
	traceback_jump origin = {{0, (unsigned int) this->sequence.size() - 1}, W_MATRIX};
	traceback_jump action;
	
	// only initize traceback if it provides free energy
	this->traceback_push(origin);///@todo use size_t
	while(this->traceback_pop(origin.target.first, origin.target.second, origin.target_matrix))
	{
#if DEBUG
		if(origin.target.first >= origin.target.second)
		{
			throw std::invalid_argument("Zuker::traceback(): invalid jump (i:" + std::to_string(origin.target.first) + " >= j:" + std::to_string(origin.target.second) + ")");
		}
#endif //DEBUG
		
		switch(origin.target_matrix)
		{
			case V_MATRIX:
				action = this->tij_v.get(origin.target);
				break;
			case W_MATRIX:
				action = this->tij_w.get(origin.target);
				break;
			case WM_MATRIX:
				action = this->tij_wm.get(origin.target);
				break;
#if DEBUG
			default:
				throw std::invalid_argument("Zuker::traceback(): (" + std::to_string(origin.target.first) + ", " + std::to_string(origin.target.second) + ") targetting from unset location\n");
				return void();
				break;
#endif //DEBUG
		}
		
		// Checks whether the action for the current position is to STORE
		if(origin.target_matrix == V_MATRIX)
		{
			this->dot_bracket.store(origin.target);
			
			independent_segment_traceback = this->sij.get(origin.target);///@todo implement it as independent_segment_traceback = this->nij.search(p); or sth like that
			if(independent_segment_traceback != nullptr)// If a Segment's traceback is found, trace its internal structure back
			{
				this->folded_segments++;
				pair_tmp = origin.target;
				
				while(independent_segment_traceback->traceback(pair_tmp.first, pair_tmp.second))
				{
#if DEBUG
					if(pair_tmp.first >= pair_tmp.second)
					{
						throw std::invalid_argument("Zuker::traceback(): segment/segmentloop introduced invalid jump (" + std::to_string(pair_tmp.first) + "," + std::to_string(pair_tmp.second) + ")\n");
					}
#endif //DEBUG
					this->dot_bracket.store(pair_tmp);
				}
			}
		}
		
		
		// For the next action (target), check whether this is a bifurcation, a threads end, or a jump to another pair
		if(action.target.first != (unsigned int) UNBOUND)
		{
			if(action.target.first == action.target.second)
			{
				if(origin.target_matrix == V_MATRIX)
				{
					this->traceback_push({{origin.target.first + 1, action.target.first}, action.target_matrix});
					this->traceback_push({{action.target.first + 1, origin.target.second - 1}, action.target_matrix});
				}
				else
				{
					this->traceback_push({{origin.target.first, action.target.first}, action.target_matrix});
					this->traceback_push({{action.target.first + 1, origin.target.second}, action.target_matrix});
				}
			}
			else
			{
				this->traceback_push(action);
			}
		}
		
		// else -> end of line / end of loop
	}
}



/**
 * @brief Pushes (i,j) & matrix-flag onto the stack
 *
 * @param arg_jump Pair to jump to + matrix to traceback from - v, w, or wm
 *
 * @todo inline
 */
void Zuker::traceback_push(traceback_jump arg_jump)
{
	this->traceback_stack.push_back(arg_jump);
}



/**
 * @brief Pops from the stack.
 *
 * @param i Nucleotide position of the sequence, paired to j, where i < j
 * @param j Nucleotide position of the sequence, paired to i, where i < j
 * @param char matrix to traceback from - v, w, or wm
 *
 * @todo use traceback_jump
 *
 * @return True for success; False otherwise.
 */
bool Zuker::traceback_pop(unsigned int &i, unsigned int &j, char &matrix)
{
	if(!this->traceback_stack.empty())
	{
		traceback_jump jump = this->traceback_stack.back();
		
		i = jump.target.first;
		j = jump.target.second;
		matrix = jump.target_matrix;
		
		this->traceback_stack.pop_back();
		
		return true;
	}
	
	return false;
}



/**
 * @brief Prints the 2D structure as DotBracket (dbn) format
 *
 * @todo Change this to Zuker::output(), add enum for OutputType::DotBracket / OutputType::ConnectivityTable / OutputType::RNAXml
 */
void Zuker::print_2D_structure(void)
{
	size_t n = this->sequence.size();
	Pair pair = Pair(0, n - 1);
	
	std::string dotbracket = "";
	this->dot_bracket.format((unsigned int) n, dotbracket); ///@todo use size_t
	
	printf(">Sequence length: %zubp, dE: %.2f kcal/mole, segments: %i\n", n, this->wij.get(pair), this->folded_segments++);
	std::cout << this->sequence.str() << "\n" << dotbracket << "\n";
}



#if DEBUG
void Zuker::_print_sij(unsigned int matrix_length)
{
	unsigned int i, j;
	
	for(i = 0; i < matrix_length; i++)
	{
		for(j = 0; j < matrix_length; j++)
		{
			Pair p = Pair(i, j);
			
			if(i > j)
			{
				std::cout << " ----    ";
			}
			else if(this->sij.get(p) == nullptr)
			{
				std::cout << " null    ";
			}
			else
			{
				std::cout << " " << this->sij.get(p);
			}
		}
		std::cout << "\n";
	}
	
	for(i = 0; i  < matrix_length * 4.4; i++)
	{
		std::cout << "--";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG



#if DEBUG
void Zuker::_print_vij(unsigned int matrix_length)
{
	unsigned int i, j;
	float p;
	
	for(i = 0; i < matrix_length; i++)
	{
		for(j = 0; j < matrix_length; j++)
		{
			Pair pair = Pair(i, j);
			
			if(i > j)
			{
				std::cout << "-------- ";
			}
			else
			{
				p = this->vij.get(pair);
				printf("%8.1f ", p);
			}
		}
		std::cout << "\n";
	}
	
	for(i = 0; i  < (matrix_length + 1) * 2; i++)
	{
		std::cout << "----";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG



#if DEBUG
void Zuker::_print_wij(unsigned int matrix_length)
{
	unsigned int i, j;
	float p;
	
	for(i = 0; i < matrix_length; i++)
	{
		for(j = 0; j < matrix_length; j++)
		{
			Pair pair = Pair(i, j);
			
			if(i > j)
			{
				std::cout << "-------- ";
			}
			else
			{
				p = this->wij.get(pair);
				printf("%8.1f ", p);
			}
		}
		std::cout << "\n";
	}
	
	for(i = 0; i  < (matrix_length + 1) * 2; i++)
	{
		std::cout << "----";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG
