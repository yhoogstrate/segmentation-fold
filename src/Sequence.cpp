/**
 * @file src/Sequence.cpp
 *
 * @date 15-apr-2015
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



#include <cstdlib>

#include "main.hpp"
#include "Nucleotide.hpp"
#include "Position.hpp"
#include "Sequence.hpp"



/**
 * @date 10-mar-2014
 */
Sequence::Sequence()
{
}



/**
 * @date 15-apr-2015
 *
 * @brief Initiates a sequences and adds all nucleoides from arg_asequence
 */
Sequence::Sequence(const char *arg_sequence)
{
	size_t i;
	
	for(i = 0; i < strlen(arg_sequence); i++)
	{
		this->push_back(arg_sequence[i]);
	}
}



/**
 * @date 18-may-2014
 *
 * @todo Make a deep-copying function if possible and write test case
 */
Sequence::Sequence(std::vector<Nucleotide> &arg_nucleotides)
{
	for(Position it = arg_nucleotides.begin(); it != arg_nucleotides.end(); ++it)
	{
		this->push_back(*it);
	}
}



/**
 * @date 18-may-2014
 */
Sequence::Sequence(std::string &arg_nucleotides)
{
	for(std::string::iterator it = arg_nucleotides.begin(); it != arg_nucleotides.end(); ++it)
	{
		this->push_back(*it);
	}
}



/**
 * @date 18-may-2014
 */
Sequence::Sequence(const std::string &arg_nucleotides)
{
	for(std::string::const_iterator it = arg_nucleotides.begin(); it != arg_nucleotides.end(); ++it)
	{
		this->push_back(*it);
	}
}



/**
 * @date 18-may-2014
 */
void Sequence::push_back(Nucleotide n)
{
	this->data.push_back(n);
}



/**
 * @date 18-may-2014
 */
void Sequence::push_back(char arg_char)
{
	switch (arg_char)
	{
		case 'a':
		case 'A':
			this->data.push_back(Nucleotide::A);
			break;
		case 'c':
		case 'C':
			this->data.push_back(Nucleotide::C);
			break;
		case 't':
		case 'T':
		case 'u':
		case 'U':
			this->data.push_back(Nucleotide::U);
			break;
		case 'g':
		case 'G':
			this->data.push_back(Nucleotide::G);
			break;
		default:
			std::string err = std::string( "Error: invalid char in sequence: '");
			err +=  arg_char;
			err +=  "' (chr " + std::to_string(arg_char) +  ")\n";
			throw std::invalid_argument( err.c_str() );
			break;
	}
}



/**
 * @brief Creates a subsequence of the sequence
 *
 * @param arg_start is the (0-based offset) base in the sequence where the subsequence starts
 * @param arg_end is the (0-based offset) base in the sequence where the subsequence ends
 *
 * @section DESCRIPTION
 * The arg_stop is more appropriate than arg_length, since the entire Zuker algorithm works with positions rather than lengths.
 *
 * @date 15-apr-2015
 */
Sequence Sequence::subseq(size_t arg_start, size_t arg_stop)
{
	std::vector<Nucleotide> vec = std::vector<Nucleotide>(this->data.begin() + arg_start, this->data.begin() + arg_stop + 1);
	return Sequence(vec);
}



/**
 * @date 18-may-2014
 * 
 * @todo check if this can be inlined?
 */
size_t Sequence::size()
{
	return this->data.size();
}



/**
 * @date 18-may-2014
 * 
 * @todo check if this can be inlined?
 */
bool Sequence::empty()
{
	return this->data.empty();
}

/**
 * @brief Obtains a requested Nucleotide within the sequence
 *
 * @section DESCRIPTION
 * The function doesn't give an out of bound error message because it
 * introduces an unnecessairy amount of comparisons; the size of the
 * sequence can be found prior to a loop using Sequence::size()
 *
 * @date 18-may-2014
 *
 * @todo templating for size_type, similar to std::vecotr::operator[] (http://www.cplusplus.com/reference/vector/vector/operator[]/)
 */
Nucleotide Sequence::operator[](size_t arg_position)
{
	return this->data[arg_position];
}



/**
 * @date 15-apr-2015
 */
std::string Sequence::str()
{
	std::string s = std::string();
	
	for(Position it = this->data.begin(); it != this->data.end(); ++it)
	{
		switch (*it)
		{
			case Nucleotide::A:
				s.push_back('A');
				break;
			case Nucleotide::C:
				s.push_back('C');
				break;
			case Nucleotide::U:
				s.push_back('U');
				break;
			case Nucleotide::G:
				s.push_back('G');
				break;
		}
	}
	
	return s;
}
