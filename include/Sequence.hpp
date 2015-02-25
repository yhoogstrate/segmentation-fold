/**
 * @file include/Sequence.hpp
 *
 * @date 25-feb-2015
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2014 Youri Hoogstrate
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



/**
 * @date 25-feb-2015
 */
class Sequence
{
	private:
		std::vector<Nucleotide> data;
		
	public:
		Sequence();
		Sequence(const char *);
		Sequence(std::string &arg_nucleotides);
		Sequence(const std::string &arg_nucleotides);
		Sequence(std::vector<Nucleotide> &arg_nucleotides);
		
		void push_back(char arg_char);
		void push_back(Nucleotide n);
		
		size_t size();
		bool empty();
		
		Sequence subseq(unsigned int arg_start, unsigned int arg_stop);
		
		std::string str();
		
		// Operators
		inline bool operator<(const Sequence &arg_left_sequence) const { return data < arg_left_sequence.data; };
		// inline bool operator>=(const Sequence &arg_left_sequence) const { return data >= arg_left_sequence.data; };
		inline bool operator==(const Sequence &arg_left_sequence) const {return data == arg_left_sequence.data;};
		Nucleotide operator[] (unsigned int);
		
		//Sequence(const Sequence& orig);
		//virtual ~Sequence();
	private:
	
};

#endif																	// SEQUENCE_HPP
