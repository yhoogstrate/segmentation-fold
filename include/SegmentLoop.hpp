/**
 * @file include/SegmentLoop.hpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2016 Youri Hoogstrate
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


#ifndef SEGMENTLOOP_HPP
#define SEGMENTLOOP_HPP


#include <array>
#include "SegmentTraceback.hpp"


/**
 * @brief Loop element that allows traceback of non-canonical pairs
 *
 * @section DESCRIPTION
 * Classical loops like triloops or tetraloops don't support the
 * possibility to include non-canonical bonds.
 *
 * <PRE>
 * 5') GGGA
 *     ||| A
 * 3') CCCA
 * </PRE>
 *
 * Imagine a loop that contains non-canonical pairs:
 *
 * <PRE>
 * 5') GGGCA UA
 *     |||:: : A
 * 3') CCCCACUA
 * </PRE>
 */
class SegmentLoop
{
	public:
		SegmentLoop(std::string arg_name, Sequence arg_sequence, std::vector<Pair> arg_bonds, float arg_gibbs_free_energy);
		
		std::string name;
		Sequence sequence;
		SegmentTraceback traceback;
		float gibbs_free_energy;
		
		Nucleotide get_nucleotide(unsigned int &arg_i);
		
		size_t size(void);
		Sequence *get_sequence(void);
};

#endif	// SEGMENTLOOP_HPP
