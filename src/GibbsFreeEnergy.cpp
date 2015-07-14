/**
 * @file src/GibbsFreeEnergy.cpp
 *
 * @date 2015-06-30
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2015 Youri Hoogstrate
 *
 * This file is part of segmentation-fold.
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


#include "main.hpp"

#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Pair.hpp"
#include "Region.hpp"

#include "Segment.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"

#include "ReadData.hpp"
#include "ReadSegments.hpp"

#include "GibbsFreeEnergy.hpp"



/**
 * @brief Constructs the GibssFreeEnergy class by interpoliting based based on the sequence size
 *
 * @date 2015-06-30
 */
GibbsFreeEnergy::GibbsFreeEnergy(Sequence &arg_sequence, ReadData &arg_thermodynamics) :
	sequence(arg_sequence),
	thermodynamics(arg_thermodynamics)
{
	this->interpolate_loop_hairpin();
	this->interpolate_loop_bulge();
	this->interpolate_loop_interior();
	
	this->interpolate_loop_hairpin_C_penalty();
}



/**
 * @brief Returns the Gibbs Free Energy of a stacking element that encloses an int11 (interior loop 1,1)
 *
 * @section DESCRIPTION
 * <PRE>
 * Pairing(i,j) Pairing(i',j') Nucleotide(i+1) Nucleotide(j-1)
 * </PRE>
 *
 * @date 2015-07-04
 */
inline float GibbsFreeEnergy::get_int11(PairingType &arg_p1, PairingType &arg_p2, Nucleotide arg_i1, Nucleotide arg_j1)
{
	return this->thermodynamics.int11[arg_p1][arg_p2][arg_i1][arg_j1];
}



/**
 * @brief Returns the Gibbs Free Energy of a stacking element that encloses an int21 (interior loop 2,1)
 *
 * @section DESCRIPTION
 * Forward:
 * <PRE>
 * Pairing(i,j)   Pairing(i',j') Nucleotide(i+1)  Nucleotide(j-1)  Nucleotide(j-2)
 * </PRE>
 *
 * Reverse:
 * <PRE>
 * Pairing(j',i') Pairing(j,  i) Nucleotide(j'+1) Nucleotide(i'-1) Nucleotide(i'-2)
 * </PRE>
 *
 * <PRE>

 5   i  i+1   i'

 3   j  j-1  j-2  j'


 * </PRE>
 *
 * @date 2015-07-04
 */
inline float GibbsFreeEnergy::get_int21(PairingType &arg_p1, PairingType &arg_p2, Nucleotide arg_i1, Nucleotide arg_j1, Nucleotide arg_j2)
{
	return this->thermodynamics.int21[arg_p1][arg_p2][arg_i1][arg_j1][arg_j2];
}



/**
 * @brief Returns the Gibbs Free Energy of a stacking element that encloses an int22
 *
 * @section DESCRIPTION
 * <PRE>
 * Pairing(i,j)   Pairing(i',j') Nucleotide(i+1)  Nucleotide(j-1)  Nucleotide(j-2)
 * Pairing(j',i') Pairing(j,  i) Nucleotide(j'+1) Nucleotide(i'-1) Nucleotide(i'-2)
 * </PRE>
 *
 * @date 2015-07-04
 */
inline float GibbsFreeEnergy::get_int22(PairingType &arg_p1, PairingType &arg_p2, Nucleotide arg_x1, Nucleotide arg_y1, Nucleotide arg_x2, Nucleotide arg_y2)
{
	return this->thermodynamics.int22[arg_p1][arg_p2][arg_x1][arg_y1][arg_x2][arg_y2];
}



/**
 * @brief Returns the Gibbs Free Energy of some predifined variables
 *
 * @date 2015-07-04
 *
 * @todo #if DEBUG -> check for out of bound
 */
inline float GibbsFreeEnergy::get_miscloop(unsigned char arg_param)
{
	return this->thermodynamics.miscloop[arg_param];
}



/**
 * @brief Returns the Gibbs Free Energy of a tloop
 *
 * @date 2015-07-04
 */
inline float GibbsFreeEnergy::get_tloop(Sequence &arg_subsequence)
{
	std::map<Sequence, float>::const_iterator search = this->thermodynamics.tloop_map.find(arg_subsequence);
	
	///@todo write as one liner () t ? f :
	if(search != this->thermodynamics.tloop_map.end())
	{
		return search->second;
	}
	else
	{
		return 0.0;
	}
}



/**
 * @brief Returns the Gibbs Free Energy of a triloop
 *
 * @date 2015-07-04
 */
inline float GibbsFreeEnergy::get_triloop(Sequence &arg_subsequence)
{
	std::map<Sequence, float>::const_iterator search = this->thermodynamics.triloop_map.find(arg_subsequence);
	
	///@todo write as oneliner () t ? f :
	if(search != this->thermodynamics.triloop_map.end())
	{
		return search->second;
	}
	else
	{
		return 0.0;
	}
}



/**
 * @brief Reterns the Gibbs Free Energy of a terminal stack of a hairpin
 *
 * @section DESCRIPTION
 * Nucleotide(i) Nucleotide(j) Nucleotide(i+1) Nucleotide(j-1)
 */
inline float GibbsFreeEnergy::get_tstackh(PairingType &arg_p1, Nucleotide arg_ip, Nucleotide arg_jp)
{
	// if debug, check for out of bound
	return this->thermodynamics.tstackh[arg_p1][arg_ip][arg_jp];
}



/**
 * @brief Returns Gibbs free energy for hairpin loops
 *
 * @date 2015-07-10
 */
inline float GibbsFreeEnergy::get_loop_hairpin(unsigned int arg_n_unpaired)
{
#if DEBUG
	///@todo Throw exception
	if(arg_n_unpaired >= (int) this->thermodynamics.loop_hairpin.size())
	{
		fprintf(stderr, "GibbsFreeEnergy::get_loop_hairpin - Out of bound: %i ,size = %i" "\n", arg_n_unpaired, (int)this->thermodynamics.loop_hairpin.size());
		exit(1);
	}
#endif //DEBUG
	
	return this->thermodynamics.loop_hairpin[arg_n_unpaired];
}



/**
 * @brief Returns Gibbs free energy for bulge loop
 *
 * @date 2015-07-10
 */
inline float GibbsFreeEnergy::get_loop_bulge(unsigned int arg_n_unpaired)
{
#if DEBUG
	///@todo Throw exception
	if(arg_n_unpaired >= (int) this->thermodynamics.loop_bulge.size())
	{
		fprintf(stderr, "GibbsFreeEnergy::get_loop_bulge - Out of bound: %i (size = %i)\n", arg_n_unpaired , (int) this->thermodynamics.loop_bulge.size());
		exit(1);
	}
#endif //DEBUG
	
	return this->thermodynamics.loop_bulge[arg_n_unpaired];
}



/**
 * @brief Returns Gibbs free energy for interior loop
 *
 * @date 2015-07-10
 */
inline float GibbsFreeEnergy::get_loop_interior(unsigned int arg_n_unpaired)
{
#if DEBUG
	///@todo Throw exception
	if(arg_n_unpaired >= (int) this->thermodynamics.loop_interior.size())
	{
		fprintf(stderr, "GibbsFreeEnergy::get_loop_interior - Out of bound: %i (size = %i)\n", arg_n_unpaired , (int) this->thermodynamics.loop_interior.size());
		exit(1);
	}
#endif //DEBUG
	
	return this->thermodynamics.loop_interior[arg_n_unpaired];
}



/**
 * @brief Returns Gibbs free energy for a poly-C loop
 *
 * @date 2015-07-10
 */
inline float GibbsFreeEnergy::get_poly_C_loop_penalty(Pair &arg_p, unsigned int n_unpaired)
{
	bool continue_loop = true;
	bool is_polyC = true;
	
	unsigned int i = arg_p.first + 1;
	unsigned int j = arg_p.second - 1;
	
	while(continue_loop)
	{
		if(this->sequence[i] != Nucleotide::C)
		{
			is_polyC = false;
		}
		
		continue_loop = (is_polyC and (i++ >= j));
	}
	
	if(is_polyC)
	{
		// if debug, check for out of bound
		return this->thermodynamics.loop_hairpin_C_penalty[n_unpaired];
	}
	else
	{
		return 0.0;
	}
}



/**
 * @brief Returns Gibbs free energy for the GGG...U loop
 *
 * @date 2015-07-10
 *
 * @todo merge this with the poly-C loop ; you can gain performance here since the one excludes the other
 */
inline float GibbsFreeEnergy::get_GGG_U_loop_penalty(Pair &arg_p)
{
	if(
		(arg_p.first >= 2) &&
		this->sequence[arg_p.first] == Nucleotide::G &&
		this->sequence[arg_p.first - 1] == Nucleotide::G &&
		this->sequence[arg_p.first - 2] == Nucleotide::G &&
		this->sequence[arg_p.second] == Nucleotide::U
	)
	{
		return this->get_miscloop(MISCLOOP_GGG_U_PENALTY);
	}
	else
	{
		return 0.0;
	}
}



/**
 * @brief Returns Gibbs free energy for connecting Stacking elements
 *
 * @date 2015-07-10
 */
inline float GibbsFreeEnergy::get_stack(Pairing &pairing1, Pairing &pairing2)
{
	return this->thermodynamics.stack[pairing1.type][pairing2.type];
}



/**
 * @brief Loads terminal stack of an interior loop
 *
 * @date 2015-06-23
 */
inline float GibbsFreeEnergy::get_tstacki(Pairing &pairing1, Nucleotide arg_i1, Nucleotide arg_j1)
{
	///@bug crashes somehow with uninitialized values
	return this->thermodynamics.tstki[pairing1.type][arg_i1][arg_j1];
}



/**
 * @brief Returns Gibbs free energy for poppen
 *
 * @date 2015-07-10
 *
 * @todo /s/poppen_p/poppen
 */
inline float GibbsFreeEnergy::get_poppen(unsigned int arg_n_unpaired_in_smallest_bulge)
{
	return this->thermodynamics.poppen_p[arg_n_unpaired_in_smallest_bulge];
}


/**
 * @brief Pre-calculates and caches all the loop_hairpin energy values based on the maximal possible length
 *
 * @date 2014-04-04
 */
inline void GibbsFreeEnergy::interpolate_loop_hairpin()
{
	unsigned int n = this->sequence.size();
	
	if(n > 30 + 2)														// if you set n = size - 2, it may become negative which is tricky since its a unsigned int
	{
		n -= 2;
		this->thermodynamics.loop_hairpin.reserve(n);					// reserve should be used, and during the for loop it should be filled using push_backs
		
		unsigned int i;
		for(i = this->thermodynamics.loop_hairpin.size(); i <= n; i++)
		{
			this->thermodynamics.loop_hairpin.push_back(this->get_loop_hairpin(30) + (this->get_miscloop(MISCLOOP_PRELOG) * log((float) i / 30)));
		}
	}
}



/**
 * @brief Caches the bulge loop element enegery values per loop size
 *
 * @section DESCRIPTION
 *
 * The largest possible bulge loop relative to the entire sequence length:
 *
 * <PRE>
 *     gcaca
 *     |  /
 * 5') G GG
 *     | | G
 * 3') C CC
 * </PRE>
 *
 * Is thus: sequence size - max hairpin size - 4
 *
 * For the given example: 12-3-4
 *
 * @date 2015-07-10
 */
inline void GibbsFreeEnergy::interpolate_loop_bulge()
{
	if(this->sequence.size() > 30 + 4 + this->thermodynamics.minimal_hairpin_length)///@todo move this comparison into the for-loop comparison -> init the loop with i = 30 ?
	{
		unsigned int i, n;
		n = this->sequence.size() - 5 - this->thermodynamics.minimal_hairpin_length;
		
		this->thermodynamics.loop_bulge.reserve(n);						// reserve should be used, and during the for loop it should be filled using push_backs
		
		for(i = this->thermodynamics.loop_bulge.size(); i <= n; i++)
		{
			this->thermodynamics.loop_bulge.push_back(this->get_loop_bulge(30) + (this->get_miscloop(MISCLOOP_PRELOG) * log((float) i / 30)));
		}
	}
}


/**
 * @brief Caches the interior loop element enegery values per loop size
 *
 * @section DESCRIPTION
 * The largest possible interior loop relative to the entire sequence length:
 *
 * <PRE>
 *     GCACA
 *     |  /
 * 5') G GG
 *     | | G
 * 3') C CC
 *     |  \
 *     CGGAG
 * </PRE>
 *
 * n_unpaired will never be longer than:
 *
 * sequence length - 4(pairs of loop) - min_hairpin_length
 *
 * @date 2015-04-30
 */
inline void GibbsFreeEnergy::interpolate_loop_interior()
{
	unsigned int accepted_diff = 5;
	
	// the following code ensures that i and n will always be >= 0 and so unsigned numbers can be used
	unsigned int i = std::max((unsigned int) 30 + 1, (unsigned int) this->thermodynamics.loop_interior.size()) + accepted_diff;
	unsigned int n = this->sequence.size() - this->thermodynamics.minimal_hairpin_length;
	
	if(n > accepted_diff)
	{
		i -= accepted_diff;
		n -= accepted_diff;
		
		this->thermodynamics.loop_interior.reserve(n);					// reserve should be used, and during the for loop it
		
		for(; i <= n; i++)
		{
			this->thermodynamics.loop_interior.push_back(
				this->thermodynamics.loop_interior[30] + (this->get_miscloop(MISCLOOP_PRELOG) * log((float) i / 30))
			);
		}
	}
}


/**
 * @brief Calculates all possible C penalty values rather than re-calculating them
 *
 * vector should be empty at startup
 * the third position (4th pos in vector) probably differs so much from interpolation that a separate variable was chosen
 *
 * @todo optimize possibly overwritten values
 *
 * @date 2015-04-04
 */
inline void GibbsFreeEnergy::interpolate_loop_hairpin_C_penalty()
{
	unsigned int n = this->sequence.size();
	unsigned int i;
	
	if(n > 2)
	{
		n -= 2;
		
		this->thermodynamics.loop_hairpin_C_penalty.reserve(n);				// reserve should be used, and during the for loop it should be filled using push_backs
		
		for(i = this->thermodynamics.loop_hairpin_C_penalty.size(); i <= n; i++)
		{
			this->thermodynamics.loop_hairpin_C_penalty.push_back(this->get_miscloop(MISCLOOP_C_HAIRPIN_INTERCEPT) + (i * this->get_miscloop(MISCLOOP_C_HAIRPIN_SLOPE)));
		}
		
		if(n >= 4)
		{
			this->thermodynamics.loop_hairpin_C_penalty[3] = this->get_miscloop(MISCLOOP_C_HAIRPIN_OF_3);
		}
	}
}




/**
 * @brief Finds the Gibbs free energy for a stacking structure element.
 *
 * <PRE>
 * 5') [i] [i+1] ...
 *      |    |
 * 3') [j] [j-1] ...
 * </PRE>
 *
 * @date 2015-02-17
 *
 * @param arg_pair Pair containing nucleotides i and j indicating their positions in the sequence, where i < j.
 *
 * @todo 3 is hardcoded? << should be used from thermodynamics?
 */
float GibbsFreeEnergy::get_stacking_pair_element(Pair &arg_pair)
{
	Nucleotide n1 = this->sequence[arg_pair.first];
	Nucleotide n2 = this->sequence[arg_pair.second];
	
	Pairing p = Pairing(n1, n2);
	
	float energy;
	
	if(((arg_pair.second - arg_pair.first - 1) < 3) || !p.is_canonical())///@todo figure out why 3 is not implemented as a constant ? shoudl be this->settings->min_hairpin_length?
	
	{
		energy = N_INFINITY;
	}
	else
	{
		Nucleotide n1p = this->sequence[arg_pair.first + 1];
		Nucleotide n2p = this->sequence[arg_pair.second - 1];
		
		Pairing p2 = Pairing(n1p, n2p);
		
		energy = this->get_stack(p, p2);///@todo figure out why in the original code this penalty was extended with eparam[0] which seemed to be a stack specific penalty of 0
	}
	
	return energy;
}


/**
 * @brief Finds the Gibbs free energy for a hairpin loop structure element.
 *
 * @section DESCRIPTION
 * The hairpinloops Gibbs free energy levels are from a fixed table upto
 * a length of 30 unpaired nucleotides. Extrapoliation is applied for a
 * number of unpaired nucleotides which >= 30.
 *
 * <PRE>
 * 5') [i] [i+1] ...
 * 3') [j] [j-1] ...
 * </PRE>
 *
 * @date 2015-02-17
 *
 * @param arg_pair Pair containing nucleotides i and j indicating their positions in the sequence, where i < j.
 *
 * @return Amount of Gibbs free energy if (i,j) enclose a hairpin loop.
 */
float GibbsFreeEnergy::get_hairpin_loop_element(Pair &arg_pair)
{
	unsigned int n_unpaired = arg_pair.second - arg_pair.first - 1;		// The number of unpaired nucleotides between i and j: if i = j+1, unPaired = 0
	float energy;
	
	if(n_unpaired < this->thermodynamics.minimal_hairpin_length)		// Too small loop
	{
		energy = N_INFINITY;///@todo this should never happen and should be taken care of in the loops
	}
	else
	{
		// Loop size + Terminal mismatch(i',j')
		//energy = this->get_loop_hairpin(n_unpaired);
		
		Nucleotide i = this->sequence[arg_pair.first];
		Nucleotide j = this->sequence[arg_pair.second];
		
		Pairing pairing = Pairing(i, j);
		
		energy = this->get_loop_hairpin(n_unpaired) +
				 this->get_tstackh(pairing.type, this->sequence[arg_pair.first + 1], this->sequence[arg_pair.second - 1]);
				 
		// Loop-sequence specific penalty
		switch(n_unpaired)
		{
			case 3:
				{
					Sequence sub_sequence = this->sequence.subseq(arg_pair.first, arg_pair.second);
					///@todo AU penalty?
					energy += this->get_triloop(sub_sequence);
					break;
				}
			case 4:
				{
					Sequence sub_sequence = this->sequence.subseq(arg_pair.first, arg_pair.second);
					energy += this->get_tloop(sub_sequence);
					break;
				}
		}
		
		energy += this->get_poly_C_loop_penalty(arg_pair, n_unpaired);
		energy += this->get_GGG_U_loop_penalty(arg_pair);
	}
	
	return energy;
}



/**
 * @brief Energy function for Bulge Loop
 *
 * @section DESCRIPTION
 * This function returns the amount of energy for a bulge loop element. The
 * parameters must be: i<ip<jp<j where in particular:
 *
 * (ip-i == 1) || (j-jp == 1)
 *
 * There are two types of bulge loops:
 *
 * - Casus "5 prime bulge":
 * <PRE>
 * 5': ... i , ... , ip ...
 * 3': ... j    ,    jp ...
 * </PRE>
 *
 * - Casus "3 prime bulge":
 * <PRE>
 * 5': ... i    ,    ip ...
 * 3': ... j , ... , jp ...
 * </PRE>
 *
 * The corresponding energy is found in additional tables.
 *
 * @date 2015-06-30
 *
 * @param i Nucleotide position of the sequence, paired to j, where i < i' (ip) < j' (jp) < j
 * @param j Nucleotide position of the sequence, paired to i, where i < i' (ip) < j' (jp) < j
 * @param ip (i') Nucleotide position of the sequence, paired to j', where i < i' (ip) < j' (jp) < j
 * @param jp (j') Nucleotide position of the sequence, paired to i', where i < i' (ip) < j' (jp) < j
 *
 * @return Amount of corresponding Gibbs free energy if there is an bulge loop betweein (i,j) and (i',j'); infinity otherwise.
 *
 *
 * @todo check: why check the size of only one size? -> does the algorithm ensure one side is always == 0?
 */
float GibbsFreeEnergy::get_bulge_loop_element(Region &r)
{
	float energy;
	
	unsigned int n_unpaired_left = r.pair2.first - r.pair1.first;
	unsigned n_unpaired_right = r.pair1.second - r.pair2.second;
	
	unsigned int n_unpaired = std::max(n_unpaired_left, n_unpaired_right) - 1;
	
	if(r.pair2.second - r.pair2.first - 1 <=  this->thermodynamics.minimal_hairpin_length)
	{
		energy = N_INFINITY;
	}
#if DEBUG
	else if(n_unpaired < 1)
	{
		printf("Warning: unneccesairy call of get_bulge_loop_element()");
		energy = N_INFINITY;
	}
#endif // DEBUG
	else
	{
		energy = this->get_loop_bulge(n_unpaired);
		
		Nucleotide n11 = this->sequence[r.pair1.first];
		Nucleotide n12 = this->sequence[r.pair1.second];
		Nucleotide n21 = this->sequence[r.pair2.first];
		Nucleotide n22 = this->sequence[r.pair2.second];
		
		Pairing pairing1 = Pairing(n11, n12);///@todo don't you get this penalty already when you calculate the hairpin?
		Pairing pairing2 = Pairing(n21, n22);
		
		if(n_unpaired == 1)
		{
			energy += this->get_stack(pairing1, pairing2);///@todo should the v matrix not be filled with the stack instead...? : return N_infinity?
		}
		else
		{
			energy += this->get_AU_penalty(pairing1);///@todo don't you get this penalty already when you calculate the hairpin?
			energy += this->get_AU_penalty(pairing2);
			
			///http://unafold.math.rpi.edu/lectures/old_RNAfold/node2.html
			// "Because free energies are assigned to loops, and not to helices, there is no a priori way of knowing whether or not a stacked pair will be terminal or not. For this reason, the terminal AU penalty is built into the TSTACKH and TSTACKI tables. For bulge, multi-branch and exterior loops, the penalty is applied explicitly. In all of these cases, the penalty is formally assigned to the adjacent loop, although it really belongs to the helix"
		}
	}
	
	return energy;
}



/**
 * @brief Energy function for Interior Loop
 *
 * @section DESCRIPTION
 * <PRE>
 * 5': i ... i' ... : 3'
 * 3': j ... j' ... : 5'
 * </PRE>
 *
 * @todo Check code for:
 * - loop contains end or beginning of sequence
 * - loop internal loop < min hairpin loop?
 *
 * @date 2013-02-17
 *
 * @param arg_region as region with Nucleotides i (pair1.first), j (pair1.second), i') and j' as positions of the sequence, where i paired to j, i' to j' and where i < i' < j' < j
 *
 * @return Amount of corresponding Gibbs free energy if there is an internal/interior loop betweein (i,j) and (i',j'); infinity otherwise.
 */
float GibbsFreeEnergy::get_interior_loop_element(Region &arg_region)
{
	unsigned int n_unpaired;
	signed int n_asym;
	
	float energy = 0.0;
	
	n_unpaired = (arg_region.pair2.first - arg_region.pair1.first) + (arg_region.pair1.second - arg_region.pair2.second) - 2;// saves one CPU operation
	n_asym = (arg_region.pair2.first - arg_region.pair1.first) - (arg_region.pair1.second - arg_region.pair2.second);// saves two CPU operations
	
	
	// 5') [i] ... [...] [i'] ... (3'
	//      |             |
	// 3') [j] ... [...] [j'] ... (5'
	if(arg_region.pair2.second - arg_region.pair2.first - 1 <=  this->thermodynamics.minimal_hairpin_length)
	{
		energy = N_INFINITY;
	}
	else if(n_unpaired > 4 || abs(n_asym) > 1)
	{
		Nucleotide ni = this->sequence[arg_region.pair1.first];
		Nucleotide nj = this->sequence[arg_region.pair1.second];
		
		Nucleotide nip = this->sequence[arg_region.pair2.first];
		Nucleotide njp = this->sequence[arg_region.pair2.second];
		
		Pairing pairing1 = Pairing(ni, nj);
		Pairing pairing2 = Pairing(njp, nip);							// I and J are switched because 5' and 3' rotate after a hairpin
		
		
		energy += this->get_loop_interior(n_unpaired);
		
		// tstki example:
		//
		//         [x1]{R1}[y2]
		// 5') xxxG            Axxxh
		//        |            |    h
		// 3') xxxC            Uxxxh
		//         [x2]{R2}[y1]
		//
		//
		// {h} and {x} are nucleotides that don't matter of which h represent the hairpin.
		//
		// G-C is the first pair (5' -> 3')
		// U-A is the second pair (5' -> 3')
		//
		// For this example you want to find the stacking panalties for forward: [G-C][x1][x2]
		// And the 'reverse' stack: [U-A][y1][y2]
		//
		// This is done in the following two lines:
		
		energy += this->get_tstacki(pairing1, this->sequence[arg_region.pair1.first + 1], this->sequence[arg_region.pair1.second - 1]);
		energy += this->get_tstacki(pairing2, this->sequence[arg_region.pair2.second + 1], this->sequence[arg_region.pair2.first - 1]); // watch out for 5' <-> 3' rotation: (jp + 1) and (ip - 1)
		
		if(n_unpaired > 30)///@todo it is not logical to have this equation only above nup=30? dbl check this with other folding algorithms >> "the f(m) array (see Ninio for details)"
		{
			energy += std::min(this->get_miscloop(MISCLOOP_ASYMETRIC_INTENRAL_LOOP), ((float) abs(n_asym) * this->get_poppen(min((int) this->thermodynamics.poppen_p.size() , min(arg_region.pair2.first - arg_region.pair1.first, arg_region.pair1.second - arg_region.pair2.second)) - 1)));
		}
	}
	else
	{
		// 5') [i] [i+1] [i']
		//      |         |
		// 3') [j] [j-1] [j']
		if(n_unpaired == 2 && n_asym == 0)
		{
			Nucleotide ni = this->sequence[arg_region.pair1.first];
			Nucleotide nj = this->sequence[arg_region.pair1.second];
			
			Nucleotide nip = this->sequence[arg_region.pair2.first];
			Nucleotide njp = this->sequence[arg_region.pair2.second];
			
			Pairing pairing1 = Pairing(ni, nj);
			Pairing pairing2 = Pairing(nip, njp);						///@todo double check this order in the .dat file
			
			
#if DEBUG
			if(pairing2.type == PairingType::None)
			{
				fprintf(stderr, "GibbsFreeEnergy::get_interior_loop_element - Requesting energy for interior loop enclosing with non canonical pairing [%i,%i].\n", arg_region.pair2.first , arg_region.pair2.second);
				energy += N_INFINITY;
			}
			else
			{
				energy += this->get_int11(pairing1.type, pairing2.type, this->sequence[arg_region.pair1.first + 1], this->sequence[arg_region.pair1.second - 1]);
			}
#else
			energy += this->get_int11(pairing1.type, pairing2.type, this->sequence[arg_region.pair1.first + 1], this->sequence[arg_region.pair1.second - 1]);
#endif //DEBUG
			
		}
		
		// 5') [i] [i+1]       [i']
		//      |               |
		// 3') [j] [j-1] [j-2] [j']
		else if(n_unpaired == 3 && n_asym == -1)
		{
			Nucleotide ni = this->sequence[arg_region.pair1.first];
			Nucleotide nj = this->sequence[arg_region.pair1.second];
			
			Nucleotide nip = this->sequence[arg_region.pair2.first];
			Nucleotide njp = this->sequence[arg_region.pair2.second];
			
			Pairing pairing1 = Pairing(ni, nj);
			Pairing pairing2 = Pairing(nip, njp);
			
			
			energy += this->get_int21(
						  pairing1.type,
						  pairing2.type,
						  this->sequence[arg_region.pair1.first + 1],
						  this->sequence[arg_region.pair1.second - 1],
						  this->sequence[arg_region.pair1.second - 2]);	// checked : correct
		}
		
		// The situation we encounter is:
		// 5') [i] [i+1] [i+2] [i'] ... (3'
		//      |               |
		// 3') [j] [j-1]       [j'] ... (5'
		//
		// If we rotate it 180 degrees, we get:
		//
		// 5') [j']       [j-1] [j] (3'
		//      |                |
		// 3') [i'] [i+2] [i+1] [i] (5'
		//
		// And we can make use of the same energy vector
		else if(n_unpaired == 3 && n_asym == 1)
		{
			Nucleotide ni = this->sequence[arg_region.pair1.first];
			Nucleotide nj = this->sequence[arg_region.pair1.second];
			
			Nucleotide nip = this->sequence[arg_region.pair2.first];
			Nucleotide njp = this->sequence[arg_region.pair2.second];
			
			Pairing pairing1 = Pairing(njp, nip);						// Rotated
			Pairing pairing2 = Pairing(nj,  ni);						// Rotated
			
			energy += this->thermodynamics.int21
					  [pairing1.type]
					  [pairing2.type]
					  [this->sequence[arg_region.pair2.second + 1]]
					  [this->sequence[arg_region.pair2.first - 1]]
					  [this->sequence[arg_region.pair2.first - 2]];
		}
		
		// 5') [i] [i+1] [i+2] [i'] (3'
		//      |               |
		// 3') [j] [j-1] [j-2] [j'] (5'
		else if(n_unpaired == 4 && n_asym == 0)//@todo figure out whether n_asym must be checked - if it would be 2, it would be a bulge loop and shouldn't reach this point. Also, the first if statement contains "|| abs(n_asym) > 1" which is always true if n_unpaired == 4 and n_asym != 0
		{
			Nucleotide ni = this->sequence[arg_region.pair1.first];
			Nucleotide nj = this->sequence[arg_region.pair1.second];
			
			Nucleotide nip = this->sequence[arg_region.pair2.first];
			Nucleotide njp = this->sequence[arg_region.pair2.second];
			
			Pairing pairing1 = Pairing(ni, nj);
			Pairing pairing2 = Pairing(nip, njp);
			
			energy += this->thermodynamics.int22
					  [pairing1.type]
					  [pairing2.type]
					  [this->sequence[arg_region.pair1.first + 1]]
					  [this->sequence[arg_region.pair1.second - 1]]
					  [this->sequence[arg_region.pair1.first + 2]]
					  [this->sequence[arg_region.pair1.second - 2]];
		}
		else//@todo figure out when this should happen - probably never
		{
			energy = N_INFINITY;
		}
	}
	
	return energy;
}



/**
 * @brief Finds a possible AU Gibbs free energy penalty for a Stacking element.
 *
 * @date 2015-02-17
 *
 * @param arg_pairing containing Nucleotide at position i of the sequence, paired to nucleotide j in the sequence, where i < j.
 *
 * @return The amount of additional Gibbs free energy as penalty for having an AU.
 */
float GibbsFreeEnergy::get_AU_penalty(Pairing &arg_pairing)
{
	if(arg_pairing.type == PairingType::AU || arg_pairing.type == PairingType::UA)
	{
		return this->get_miscloop(MISCLOOP_AU_PENALTY);
	}
	else
	{
		return 0.0;
	}
}



/**
 * @brief Finds an estimated amount of Gibbs free energy for a canonical pairing between two nucleotides when surrounded by a segment.
 *
 * @todo move this to config file
 *
 * @date 2015-02-20
 *
 * @todo change into stacking penalty or something - penalties are independent from the neighbours?
 * @todo check if this data is not in one of the other stacking files
 * - otherwise represent them as functions (this->get_stack(p.first,p.second) + this->get_stack(p.second,p.first)/2) and add some loading/interpolation function to avoid many calculations
 * - OR WHETHER THIS IS NOT IDENTICAL TO thermodynamics.stack_p[pairing1,pairing2]
 */
float GibbsFreeEnergy::get_stacking_pair_without_surrounding(PairingPlus &arg_p)
{
	float energy;
	
	switch(arg_p.type)
	{
		case PairingType::GC:
			energy = -2.55;
			break;
		case PairingType::CG:
			energy = -2.10;
			break;
			
		case PairingType::UA:
			energy = -1.50;
			break;
		case PairingType::AU:
			energy = -1.38;
			break;
			
		case PairingType::GU:
			energy = -1.01;
			break;
		case PairingType::UG:
			energy = -0.78;
			break;
		default:
			energy = N_INFINITY;
			break;
	}
	
	return energy;
}
