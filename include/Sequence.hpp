/**
 * @file include/Sequence.hpp
 *
 * @date 2015-12-06
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



#ifndef SEQUENCE_HPP
#define	SEQUENCE_HPP



#include <string.h>
#include <vector>
#include <algorithm>

#include "SubSequence.hpp"


/**
 * @brief A RNA or DNA Sequence object, primarily used as the Sequence to be folded
 *
 * @date 2015-12-05
 */
class Sequence
{
	public:
	
		/// The sequence as vector of nucleotides
		std::vector<Nucleotide> data;
		
		Sequence();
		Sequence(const char *);
		Sequence(std::string &arg_nucleotides);
		Sequence(const std::string &arg_nucleotides);
		Sequence(std::vector<Nucleotide> &arg_nucleotides);
		
		void push_back(char arg_char);
		void push_back(Nucleotide n);
		
		size_t size();
		bool empty();
		
		Sequence subseq(size_t arg_start, size_t arg_stop);
		SubSequence ssubseq(size_t arg_start, size_t arg_stop);
		
		std::string str();
		
		
		/// Smaller than operator that compares whether the vector of Nucleotides is smaller than the other
		inline bool operator<(const Sequence &arg_left_sequence) const { return data < arg_left_sequence.data;};
		/// The equality operator that returns equal only if the sequences have the same composition
		inline bool operator==(const Sequence &arg_left_sequence) const {return data == arg_left_sequence.data;};
		
		
		///@note The following line is not being used in the code, but might be desired to include it as API
		// inline bool operator>=(const Sequence &arg_left_sequence) const { return data >= arg_left_sequence.data;};
		
		char compare(Sequence &arg_query);//    returns IS_EQUAL, IS_SMALLER (if this is smaller) or IS_LARGER (if arg_query is smaller)
		char compare(SubSequence &arg_query);// returns IS_EQUAL, IS_SMALLER (if this is smaller) or IS_LARGER (if arg_query is smaller)
		
		Nucleotide operator[](size_t);
};

#endif	// SEQUENCE_HPP
