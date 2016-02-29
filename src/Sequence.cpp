/**
 * @file src/Sequence.cpp
 *
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



#include "main.hpp"
#include "Nucleotide.hpp"
#include "Position.hpp"
#include "Sequence.hpp"



/**
 */
Sequence::Sequence()
{
}



/**
 * @brief Initiates a sequences and adds all nucleoides from arg_asequence
 *
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
 * @brief Initiates a sequences and adds all nucleoides from arg_asequence
 *
 */
Sequence::Sequence(std::vector<Nucleotide> &arg_nucleotides)
{
	for(Position iterator = arg_nucleotides.begin(); iterator != arg_nucleotides.end(); ++iterator)
	{
		this->push_back(*iterator);
	}
}



/**
 * @brief Initializes a Sequence from a std::string of charset "^[actguACTUG ]+$"
 *
 */
Sequence::Sequence(std::string &arg_nucleotides)
{
	for(std::string::iterator iterator = arg_nucleotides.begin(); iterator != arg_nucleotides.end(); ++iterator)
	{
		this->push_back(*iterator);
	}
}



/**
 * @brief Initializes a Sequence from a const std::string of charset "^[actguACTUG ]+$"
 *
 */
Sequence::Sequence(const std::string &arg_nucleotides)
{
	for(std::string::const_iterator iterator = arg_nucleotides.begin(); iterator != arg_nucleotides.end(); ++iterator)
	{
		this->push_back(*iterator);
	}
}



/**
 * @brief Adds a Nucleotide to the end of the Sequence
 *
 */
void Sequence::push_back(Nucleotide nucleotide)
{
	this->data.push_back(nucleotide);
}



/**
 * @brief Adds a char of charset "^[actguACTUG]$" as Nucleotide to the Sequence
 *
 */
void Sequence::push_back(char arg_char)
{
	switch(arg_char)
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
			std::string err = std::string("Error: invalid char in sequence: '");
			err +=  arg_char;
			err +=  "' (chr " + std::to_string(arg_char) +  ")\n";
			throw std::invalid_argument(err.c_str());
			break;
	}
}



/**
 * @brief Creates a subsequence of the sequence
 *
 * @section DESCRIPTION
 * The arg_stop is more appropriate than arg_length, since the entire Zuker algorithm works with positions rather than lengths.
 *
 * @param arg_start is the (0-based offset) nucleotide in the sequence where the subsequence starts
 * @param arg_stop is the (0-based offset) nucleotide in the sequence where the subsequence ends
 *
 */
Sequence Sequence::subseq(size_t arg_start, size_t arg_stop)
{
	std::vector<Nucleotide> vec = std::vector<Nucleotide>(this->data.begin() + (long) arg_start, this->data.begin() + (long) arg_stop + 1);
	return Sequence(vec);
}



/**
 * @brief Creates a subsequence of the sequence
 *
 * @section DESCRIPTION
 * The arg_stop is more appropriate than arg_length, since the entire Zuker algorithm works with positions rather than lengths.
 *
 * @param arg_start is the (0-based offset) nucleotide in the sequence where the subsequence starts
 * @param arg_stop is the (0-based offset) nucleotide in the sequence where the subsequence ends
 *
 */
SubSequence Sequence::ssubseq(size_t arg_start, size_t arg_stop)
{
	return SubSequence(this->data.begin() + (long) arg_start, this->data.begin() + (long) arg_stop);
}



/**
 * @brief Gives the number of Nucleotides in the Sequence
 *
 */
size_t Sequence::size()
{
	return this->data.size();
}



/**
 * @brief Returns whether the Sequence is empty or not
 *
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
 */
Nucleotide Sequence::operator[](size_t arg_position)
{
	return this->data[arg_position];
}



/**
 * @brief Converts the sequence to a std::string
 *
 * @section DESCRIPTION
 * Nucleotides are represented in uppercase and U's are chosen over T's
 * because it's in an RNA context
 *
 */
std::string Sequence::str()
{
	std::string sequence = std::string();
	
	for(Position it = this->data.begin(); it != this->data.end(); ++it)
	{
		switch(*it)
		{
			case Nucleotide::A:
				sequence.push_back('A');
				break;
			case Nucleotide::C:
				sequence.push_back('C');
				break;
			case Nucleotide::U:
				sequence.push_back('U');
				break;
			case Nucleotide::G:
				sequence.push_back('G');
				break;
		}
	}
	
	return sequence;
}



/**
 * @brief returns IS_EQUAL, IS_SMALLER (if this is smaller) or IS_LARGER (if arg_query is smaller)
 *
 */
char Sequence::compare(Sequence &arg_query)
{
	size_t size_this = this->size();
	size_t size_query = arg_query.size();
	
	if(size_this < size_query)
	{
		return IS_SMALLER;
	}
	else if(size_this == size_query)
	{
		for(unsigned int i = 0; i < size_this; i++)
		{
			if(this->data[i] < arg_query[i])
			{
				return IS_SMALLER;
			}
			else if(this->data[i] > arg_query[i])
			{
				return IS_LARGER;
			}
		}
		
		return IS_EQUAL;
	}
	else
	{
		return IS_LARGER;
	}
}



/**
 * @brief Compares whether this sequence is smaller, larger or equal compared to a given SubSequence
 *
 */
char Sequence::compare(SubSequence &arg_query)
{
	size_t size_this = this->size();
	
	if(size_this < arg_query.size)
	{
		return IS_SMALLER;
	}
	else if(size_this == arg_query.size)
	{
		for(size_t i = 0; i < size_this; i++)
		{
			Nucleotide n_i = arg_query[i];
			if(this->data[i] < n_i)
			{
				return IS_SMALLER;
			}
			else if(this->data[i] > n_i)
			{
				return IS_LARGER;
			}
		}
		
		return IS_EQUAL;
	}
	else
	{
		return IS_LARGER;
	}
	
	return IS_EQUAL;
}
