/**
 * @file include/GibbsFreeEnergy.hpp
 *
 * @date 2015-06-30
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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
 * </PRE>
 */

#ifndef GIBBSFREEENERGY_HPP
#define	GIBBSFREEENERGY_HPP


class GibbsFreeEnergy
{
	private:
		inline float get_int11(PairingType &arg_p1, PairingType &arg_p2, Nucleotide arg_i1, Nucleotide arg_j1);/// @note historically known as 'sint2()'
		inline float get_int21(PairingType &arg_p1, PairingType &arg_p2, Nucleotide arg_i1, Nucleotide arg_j1, Nucleotide arg_j2);/// @note historically known as  'asint1x2()' and 'asint3()'
		inline float get_int22(PairingType &arg_p1, PairingType &arg_p2, Nucleotide arg_i1, Nucleotide arg_j1, Nucleotide arg_i2, Nucleotide arg_j2);/// @note historically known as  'sint4()'
		
		inline float get_miscloop(unsigned char arg_param);
		
		/// @todo merge these three functions into get_loop and use the sizes to further distinguish between tetra and triloop
		inline float get_triloop(Sequence &arg_subsequence);			///@todo consider using &SubSequence as argument
		inline float get_tloop(Sequence &arg_subsequence);
		
		inline float get_stack(Pairing &paring1, Pairing &pairing2);// Stack after stack
		inline float get_tstacki(Pairing &paring1, Nucleotide arg_i1, Nucleotide arg_j1);// Terminal stack interior
		inline float get_tstackh(PairingType &arg_p1, Nucleotide arg_x2, Nucleotide arg_y2);// Terminal stack hairpin
		
		inline float get_loop_hairpin(unsigned int arg_n_unpaired);
		inline float get_loop_bulge(unsigned int arg_n_unpaired);
		inline float get_loop_interior(unsigned int arg_n_unpaired);
		
		inline float get_poppen(unsigned int arg_n_unpaired_in_smallest_bulge);
		
		inline float get_poly_C_loop_penalty(Pair &arg_p, unsigned int arg_n_unpaired);
		inline float get_GGG_U_loop_penalty(Pair &arg_p);
		
		inline void interpolate_loop_hairpin();
		inline void interpolate_loop_bulge();
		inline void interpolate_loop_interior();
		inline void interpolate_loop_hairpin_C_penalty();
		
		//@note the following functions might me implemented one day
		//inline float get_tstackcoax();
		//inline float get_tstackm();
		//inline float get_tstack();
		//inline float get_coaxstack();
		//inline float get_coaxial();
		//inline float get_dangle();
		
	protected:
		Sequence &sequence;
		ReadData &thermodynamics;
		
	public:
		GibbsFreeEnergy(Sequence &arg_sequence, ReadData &arg_thermodynamics);
		
		float get_stacking_pair_element(Pair &p);
		float get_hairpin_loop_element(Pair &p);
		float get_bulge_loop_element(Region &r);
		float get_interior_loop_element(Region &r);
		
		float get_AU_penalty(Pairing &p);
		float get_stacking_pair_without_surrounding(PairingPlus &p);
		
		Segment *get_segment(Sequence &arg_sequence_5p, Sequence &arg_sequence_3p);
};

#endif	// GIBBSFREEENERGY_HPP

