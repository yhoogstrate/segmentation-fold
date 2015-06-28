/**
 * @file include/SubSequence.hpp
 *
 * @date 2015-05-06
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


#ifndef SUBSEQUENCE_HPP
#define	SUBSEQUENCE_HPP


#include "main.hpp"

#include "Nucleotide.hpp"
#include "Position.hpp"


/**
 * @brief
 *
 * @date 2015-05-06
 */
class SubSequence
{
#if DEBUG
	private:
		void _check_order(void);
#endif //DEBUG
		
	public:
		Position position1;
		Position position2;
		
		size_t size;/// The length of the subsequence [18,18] << subsequence of only one nucleotide has size 1
		
		SubSequence(Position arg_position1, Position arg_position2);
		SubSequence(Position arg_position1, Position arg_position2, size_t arg_nucleotides_inbetween);
		
		void init(void);
		
		Nucleotide operator[](size_t);
};

#endif	// SUBSEQUENCE_HPP
