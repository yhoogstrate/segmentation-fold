/**
 * @file src/SegmentLoop.cpp
 *
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 * </PRE>
 */



#include "main.hpp"

#include <vector>

#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "SegmentLoop.hpp"



/**
 * @brief Constructor of SegmentLoop class.
 *
 */
SegmentLoop::SegmentLoop(std::string arg_name, Sequence arg_sequence, std::vector<Pair> arg_bonds, float arg_gibbs_free_energy):
	name(arg_name),
	sequence(arg_sequence),
	traceback(arg_bonds),
	gibbs_free_energy(arg_gibbs_free_energy)
{
}



/**
 * @brief Return the size of the segmentloops internal sequence
 *
 * */
size_t SegmentLoop::size(void)
{
	return this->sequence.size();
}



/**
 * @brief Returns a particular nucleotide in the 2 sequence of the segmentloop
 *
 * @param i Location (starting from 0, ending at n-1) in the sequence of interest
 *
 */
Nucleotide SegmentLoop::get_nucleotide(unsigned int &i)
{
	return this->sequence[i];
}



/**
 * @brief Returns sequences of the segmentloop
 *
 */
Sequence *SegmentLoop::get_sequence(void)
{
	return &this->sequence;
}
