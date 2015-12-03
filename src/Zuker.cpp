/**
 * @file src/Zuker.cpp
 *
 * @date 2015-12-03
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2015 Youri Hoogstrate
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
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"

#include "ScoringTree.hpp"

#include "ScoringMatrix.hpp"
#include "Settings.hpp"
#include "DotBracket.hpp"
#include "ReadData.hpp"

#include "Zuker.hpp"



/**
 * @brief Constructs /initializes the Zuker class: include parameters and init an empty dotbracket output.
 *
 * @date 2012-11-05
 *
 * @todo move this to this->init(); and run this->init(); or rename it to this->reset();
 */
Zuker::Zuker(Settings &arg_settings, Sequence &arg_sequence, ReadData &arg_thermodynamics) :
	GibbsFreeEnergy(arg_sequence, arg_thermodynamics),
	settings(arg_settings),
	pij(arg_sequence.size(), UNBOUND),
	qij(arg_sequence.size(), UNBOUND),
	vij(arg_sequence.size(), N_INFINITY),
	wij(arg_sequence.size(), 0.0),
	pathmatrix_corrected_from(arg_sequence.size(), false),
	loopmatrix(arg_sequence.size(), Pair(0, 0)),
	sij(arg_sequence.size(), nullptr)
{
	this->pij.fill(NOT_YET_CALCULATED);
	
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
 * @date 2013-09-20
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
 * @date 2015-12-03
 *
 * @param p1 A pair of positions refering to Nucleotide positions in the sequence, where pi.first < p1.second
 *
 * @return amount of Gibbs free energy provided for folding nucleotide i with j assuming i and j are paired
 */
float Zuker::v(Pair &p1, PairingPlus &p1p)
{
	unsigned int ip;
	unsigned int jp;
	
	float energy, tmp, tmp_k;
	
	Segment *tmp_segment = nullptr;
	Segment *tmp_segment2;
	Pair tmp_loopmatrix_value;
	
	if(this->pij.get(p1) > NOT_YET_CALCULATED)			///@todo -create bool function > {pij}.is_calculated()    @todo2 use max unsigned int value  // This point is already calculated; efficiency of dynamic programming
	{
		energy = this->vij.get(p1);
	}
	else
	{
	
#if DEBUG
		if(!p1p.is_canonical() || (p1.second - p1.first) <= this->settings.minimal_hairpin_length)
		{
			throw std::invalid_argument("Zuker::(" + std::to_string(p1.first) + ", " + std::to_string(p1.second) + "): this pair should never be checked within this function because it's energy is infinity by definition");
		}
		else
		{
#endif //DEBUG
		
			energy = this->get_hairpin_loop_element(p1);				// Hairpin element
			tmp_loopmatrix_value = {p1.first + 1, p1.second - 1};
			
			for(ip = p1.first + 1; ip < p1.second; ip++)
			{
				for(jp = p1.second - 1; jp > ip; jp--)
				{
					PairingPlus pairing2 = PairingPlus(this->sequence_begin + ip, this->sequence_begin + jp);
					if(pairing2.type != PairingType::None)				// The following structure elements must be enclosed by pairings on both sides
					{
						Pair p2 = Pair(ip, jp);
						float v_ij_jp = this->v(p2, pairing2);
						
						Region region = Region {p1, p2};
						
						if(ip == (p1.first + 1) && jp == (p1.second - 1))// Stacking element
						{
							tmp = this->get_stacking_pair_element(p1) + v_ij_jp;
							if(tmp < energy)
							{
								energy = tmp;
								
								tmp_loopmatrix_value = {p1.first + 1, p1.second - 1};
								tmp_segment = nullptr;
							}
						}
						else if(ip == (p1.first + 1) || jp == (p1.second - 1))// Bulge-loop element
						{
							tmp = this->get_bulge_loop_element(region) + v_ij_jp;
							
							if(tmp < energy)
							{
								energy = tmp;
								
								tmp_loopmatrix_value = {ip, jp};
								tmp_segment = nullptr;
							}
						}
						else											// Interior loop or Segment
						{
							tmp = this->get_interior_loop_element(region) + v_ij_jp;
							
							if(tmp < energy)
							{
								energy = tmp;
								
								tmp_loopmatrix_value = {ip, jp};
								tmp_segment = nullptr;
							}
							
							// Find segments:
							SubSequence pp1 = SubSequence(this->sequence_begin + p1.first + 1, this->sequence_begin + ip - 1);
							SubSequence pp2 = SubSequence(this->sequence_begin + jp + 1, this->sequence_begin + p1.second - 1);
							tmp_segment2 = this->thermodynamics.segments.search(pp1 , pp2);
							
							if(tmp_segment2 != nullptr)
							{
								//Nucleotide n1 = this->sequence[p1.first];
								//Nucleotide n2 = this->sequence[p1.second];
								//Pairing pairing = Pairing(n1, n2);
								
								tmp_k = tmp_segment2->gibbs_free_energy + this->get_stacking_pair_without_surrounding(p1p) + v_ij_jp;
							}
							else
							{
								tmp_k = N_INFINITY;
							}
							if(tmp_k < energy)
							{
								energy = tmp_k;
								
								tmp_loopmatrix_value = {ip, jp};
								tmp_segment = tmp_segment2;
							}
						}
					}
				}
				
				/*
				i=0, j=11 < indicated in brackets
				     *       i' = 5
				[(...)(...)]
				
				  *          i' = 2
				[()(......)]
				
				        *    i' = 8
				[(......)()]
				
				-> earlier versions compromised this function
				   by using it to 'extend' their stack
				   one bell of the fork was 0 size and the other bell
				   was the remainder
				*/
				if(ip > p1.first + 1 && ip < p1.second - 2)												// Multi-loop element
				{
					Pair p3 = Pair(p1.first + 1, ip);
					Pair p4 = Pair(ip + 1, p1.second - 1);
					Region region = Region {p3, p4};
					tmp = this->energy_bifurcation(region);//@todo energy paired bifurcation?
					
					if(tmp < energy)
					{
						energy = tmp;
						
						tmp_loopmatrix_value = {ip, jp};
						tmp_segment = nullptr;
					}
				}
			}
			
			this->loopmatrix.set(p1.first, p1.second, tmp_loopmatrix_value);
			if(tmp_segment != nullptr)
			{
				this->sij.set(p1, &tmp_segment->traceback);
			}
#if DEBUG
		}
#endif //DEBUG
		
		this->vij.set(p1.first, p1.second, energy);
	}
	
	return energy;
}



/**
 * @brief Wij Function - provided Gibbs free energy for sequence i...j
 *
 * @date 2015-12-01
 *
 * @param p1 A pair of positions refering to Nucleotide positions in the sequence, where pi.first < p1.second
 *
 * @return amount of Gibbs free energy provided for folding nucleotide i with j
 */
float Zuker::w(Pair &p1)
{
	float energy;
	int tmp_pij, tmp_qij;
	
	if(this->pij.get(p1) > NOT_YET_CALCULATED)							// Already calculated or fixed by initialization; dynamic programming function
	{
		return this->wij.get(p1);
	}
	else
	{
		int k = 0;///@todo check whether it can be unset
		float tmp;
		bool tmp_path_matrix = 0;
		
		if((p1.second - p1.first) <= (this->settings.minimal_hairpin_length))
		{
			this->vij.set(p1, N_INFINITY);
			energy = 0.0;
			
			tmp_pij = UNBOUND;
			tmp_qij = UNBOUND;
		}
		else
		{
			tmp_path_matrix = 1;
			
			PairingPlus p1p = PairingPlus(this->sequence_begin + p1.first, this->sequence_begin + p1.second);
			
			if(!p1p.is_canonical())
			{
				tmp_pij = UNBOUND;
				tmp_qij = UNBOUND;
				
				energy = N_INFINITY;
			}
			else
			{
				tmp_pij = BOUND;
				tmp_qij = BOUND;
				
				energy = this->v(p1, p1p);
			}
			
			/*
			following extreme w()-directed bifurcations are possible
			
			 *       k = i + 1
			[)(....]
			
			     *   k = j - 2; k < j - 1
			[....)(]
			 */
			for(k = p1.first + 1; k < p1.second - 1; k++)				// Find bifurcation in non-paired region
			{
				Pair p2 = Pair(p1.first, k);
				Pair p3 = Pair(k + 1, p1.second);
				
				tmp = this->w(p2) + this->w(p3);
				
				if(tmp < energy)
				{
					// Can also be done by checking whether pathmatrix_corrected_from > 0 ? >> and only store those positions in a tree instead of an entire matrix
					tmp_path_matrix = 0;
					
					energy = tmp;
					tmp_pij = k;
				}
			}
			
		}
		
		///@todo get some way to obtain the index just once? -> this could be generalized at the top level as well (add it to the pair for example)
		this->pathmatrix_corrected_from.set(p1, tmp_path_matrix);
		this->pij.set(p1, tmp_pij);
		this->qij.set(p1, tmp_qij);
		this->wij.set(p1, energy);
	}
	
	return energy;
}



/**
 * @brief Calculates the energy of a fork / bifurcation
 *
 * @section DESCRIPTION
 * This function is part of the Zuker class because its an energy
 * function using only values from the w-matrix
 *
 * @date 2015-04-03
 */
inline float Zuker::energy_bifurcation(Region &region)
{
	return this->w(region.pair1) + this->w(region.pair2);
}


/**
 * @brief The traceback algorithm, finds the optimal path through the matrices.
 *
 * @date 2015-12-01
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
	int i, j, k, ip, jp;
	int tmp_i, tmp_j;
	SegmentTraceback *independent_segment_traceback;
	
	bool pick_from_v_path;
	
	this->traceback_push(0, this->sequence.size() - 1, false);
	
	while(this->traceback_pop(&i, &j, &pick_from_v_path))
	{
#if DEBUG
		if(i < j)
		{
#endif //DEBUG
			Pair pair1 = Pair(i, j);
			k = this->pij.get(pair1);
			
			Pair pair2 = this->loopmatrix.get(pair1); // continue with Pair object instead
			ip = pair2.first;
			jp = pair2.second;
			
			// [if from the v path ] or [from a fork; V or W fork?]
			if(pick_from_v_path == true || this->pathmatrix_corrected_from.get(pair1) != 0)	// Decide which matrix to pick from
			
				/** @todo check whether it works whenever the FIRST fold is a MOTIF */
			{
				k = this->qij.get(pair1);
				pick_from_v_path = true;
			}
			
			
			if(k == BOUND)												// Parse the route
			{
				this->dot_bracket.store(i, j);
				
				independent_segment_traceback = this->sij.get(pair1);					///@todo implement it as independent_segment_traceback = this->nij.search(p); or sth like that
				
				if(independent_segment_traceback != nullptr)								// If a Segment's traceback is found, trace its internal structure back
				{
					tmp_i = i;
					tmp_j = j;
					
					while(independent_segment_traceback->traceback(tmp_i, tmp_j))
					{
#if DEBUG
						if(tmp_j <= tmp_i)
						{
							throw std::invalid_argument("Traceback with segment introduced incorrect jump (" + std::to_string(tmp_i) + "," + std::to_string(tmp_j) + ")\n");
						}
#endif //DEBUG
						this->dot_bracket.store(tmp_i, tmp_j);
					}
				}
				
				// Continue with enclosing base pair.
				ip = this->loopmatrix.get(pair1).first;
				jp = this->loopmatrix.get(pair1).second;
				
				if(ip != jp)
				{
					this->traceback_push(ip, jp, pick_from_v_path);
				}
				else													// fork from paired; v()
				{
#if DEBUG
					if(ip <= i + 1)
					{
						throw std::invalid_argument("Traceback introduced incorrect jump from paired bifurcation (i+1 == i'): " + std::to_string(i + 1) + "," + std::to_string(ip) + "\n");
					}
					
					if(ip >= j - 2)
					{
						throw std::invalid_argument("Traceback introduced incorrect jump from paired bifurcation (i'+1 == j-1): " + std::to_string(ip + 1) + "," + std::to_string(j = 1) + "\n");
					}
#endif //DEBUG
					// There are two loops within the current pair
					
					this->traceback_push(i + 1, ip, false);
					this->traceback_push(ip + 1, j - 1, false);
				}
			}
			else if(k >= 0)												// fork from non paired; w()
			{
#if DEBUG
				if(k <= i)
				{
					throw std::invalid_argument("Traceback introduced incorrect jump from unpaired bifurcation (i+1 == i'): " + std::to_string(i + 1) + "," + std::to_string(k) + "\n");
				}
				
				if(k >= j - 1)
				{
					throw std::invalid_argument("Traceback introduced incorrect jump from unpaired bifurcation (i'+1 == j-1): " + std::to_string(k + 1) + "," + std::to_string(j = 1) + "\n");
				}
#endif //DEBUG
				// The current pair actually forms 2 loops
				
				this->traceback_push(i, k, false);
				this->traceback_push(k + 1, j, false);
			}
#if DEBUG
		}
		else
		{
			throw std::invalid_argument("Traceback encountered an incorrect jump (i:" + std::to_string(i) + " >= j:" + std::to_string(j) + ")");
		}
#endif //DEBUG
	}
}



/**
 * @brief Pushes (i,j) & matrix-flag onto the stack
 *
 * @date 2012-12-06
 *
 * @param i Nucleotide position of the sequence, paired to j, where i < j
 * @param j Nucleotide position of the sequence, paired to i, where i < j
 * @param pick_from_v_path
 *
 * @todo use Pair() instead of i and j
 */
void Zuker::traceback_push(int i, int j, bool pick_from_v_path)
{
	this->traceback_stack.push_back({i, j, pick_from_v_path});
}



/**
 * @brief Pops from the stack.
 *
 * @date 2013-10-02
 *
 * @param i Nucleotide position of the sequence, paired to j, where i < j
 * @param j Nucleotide position of the sequence, paired to i, where i < j
 * @param pick_from_v_path
 *
 * @return True for success; False otherwise.
 */
bool Zuker::traceback_pop(int *i, int *j, bool *pick_from_v_path)
{
	if(not this->traceback_stack.empty())
	{
		traceback_jump jump = this->traceback_stack.back();
		
		(*i) = jump.i;
		(*j) = jump.j;
		
		(*pick_from_v_path) = jump.jump_to_v_path;
		
		this->traceback_stack.pop_back();
		
		return true;
	}
	
	return false;
}



/**
 * @brief Prints the 2D structure as DotBracket (dbn) format
 *
 * @date 2015-07-13
 *
 * @todo Change this to Zuker::output(), add enum for OutputType::DotBracket / OutputType::ConnectivityTable / OutputType::RNAXml
 */
void Zuker::print_2D_structure(void)
{
	size_t n = this->sequence.size();
	Pair pair = Pair(0, n - 1);
	
	std::string dotbracket = "";
	this->dot_bracket.format(n, dotbracket);
	
	printf(">Sequence length: %zubp, dE: %g kcal/mole\n", n, this->wij.get(pair));
	printf("%s\n", this->sequence.str().c_str());
	std::cout << dotbracket << "\n";
}



#if DEBUG
void Zuker::_print_pathmatrix_corrected_from(unsigned int matrix_length)
{
	for(int i = 0; i  < matrix_length; i++)
	{
		for(int j = 0; j  < matrix_length; j++)
		{
			Pair p = Pair(i, j);
			if(i > j)
			{
				std::cout << " -";
			}
			else if(this->pathmatrix_corrected_from.get(p) != 0)
			{
				std::cout << " t";
			}
			else
			{
				std::cout << " f";
			}
		}
		
		std::cout << "\n";
	}
	
	for(int i = 0; i < matrix_length; i++)
	{
		std::cout << "--";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG



#if DEBUG
void Zuker::_print_loopmatrix(unsigned int matrix_length)
{
	for(int i = 0; i  < matrix_length; i++)
	{
		for(int j = 0; j  < matrix_length; j++)
		{
			Pair p1 = Pair(i, j);
			if(i > j)
			{
				std::cout << " -,-";
			}
			else
			{
				Pair p = this->loopmatrix.get(p1);
				std::cout << " " << p.first << "," << p.second;
			}
		}
		std::cout << "\n";
	}
	
	for(int i = 0; i  < matrix_length * 2; i++)
	{
		std::cout << "--";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG



#if DEBUG
void Zuker::_print_nij(unsigned int matrix_length)
{
	for(int i = 0; i  < matrix_length; i++)
	{
		for(int j = 0; j  < matrix_length; j++)
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
	
	for(int i = 0; i  < matrix_length * 4.4; i++)
	{
		std::cout << "--";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG



#if DEBUG
void Zuker::_print_pij(unsigned int matrix_length)
{
	int p;
	
	for(int i = 0; i < matrix_length; i++)
	{
		for(int j = 0; j < matrix_length; j++)
		{
			Pair pair = Pair(i, j);
			
			if(i > j)
			{
				std::cout << " ---";
			}
			else
			{
				p = this->pij.get(pair);
				
				std::cout << " ";
				if(p >= 0)
				{
					std::cout << " ";
				}
				if(p > -10 && p < 10)
				{
					std::cout << " ";
				}
				std::cout << p;
			}
		}
		std::cout << "\n";
	}
	
	for(int i = 0; i  < matrix_length * 2; i++)
	{
		std::cout << "--";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG



#if DEBUG
void Zuker::_print_vij(unsigned int matrix_length)
{
	float p;
	
	for(int i = 0; i < matrix_length; i++)
	{
		for(int j = 0; j < matrix_length; j++)
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
	
	for(int i = 0; i  < (matrix_length + 1) * 2; i++)
	{
		std::cout << "----";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG



#if DEBUG
void Zuker::_print_wij(unsigned int matrix_length)
{
	float p;
	
	for(int i = 0; i < matrix_length; i++)
	{
		for(int j = 0; j < matrix_length; j++)
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
	
	for(int i = 0; i  < (matrix_length + 1) * 2; i++)
	{
		std::cout << "----";
	}
	
	std::cout << "\n\n";
}
#endif // DEBUG
