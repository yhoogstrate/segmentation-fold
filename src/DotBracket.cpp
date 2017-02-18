/**
 * @file src/DotBracket.cpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2016 Youri Hoogstrate
 *
 * This file is part of segmentation-fold and originally taken from
 * yh-kt-fold.
 *
 * segmentation-fold is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
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

#include "Utils/utils.hpp"
#include "DotBracket.hpp"


DotBracket::DotBracket()
{
}



/**
 * @brief Searches for nucleotide pairing with i
 *
 * @param i position i in the RNA sequence
 * @return The position in the sequence pairing with the nucleotide at position i, "UNBOUND" otherwise
 *
 * @todo Make this a WHILE loop and avoid the early returns
 * @todo get rid of signed here - maybe make an at() function that requires to find something
 */
signed int DotBracket::find(unsigned int arg_i)
{
	unsigned int i;
	
	for(i = 0; i < this->pairings.size(); i++)
	{
		if((unsigned int) this->pairings[i].first == arg_i)
		{
			return (signed int) this->pairings[i].second;
		}
		else if((unsigned int)  this->pairings[i].second == arg_i)
		{
			return (signed int) this->pairings[i].first;
		}
	}
	
	return UNBOUND;
}



/**
 * @brief Stores pairing (i,j)
 *
 * @param arg_pair Nucleotides position in the RNA sequence
 */
void DotBracket::store(Pair &arg_pair)
{
	this->pairings.push_back(arg_pair);
}



/**
 * @brief Formats and prints the 2D structure in DotBracket format
 *
 * @param rnaSequenceLength The length of the sequence of which the 2D strucutre is calculated
 *
 * @todo Make it work with streams, so that output can be put in a file.
 * @todo Make the variables private const static class variables
 */
void DotBracket::format(unsigned int n, std::string &output)
{
	///@todo -> size_t or template
	signed int i;
	signed int j;
	
	for(i = 0; i < (signed int) n; i++)
	{
		j = this->find((unsigned int) i);
		
		if(j == UNBOUND)
		{
			output.push_back(DOTBRACKET__NO_PAIRING);					// '.'
		}
		else if(j > i)
		{
			output.push_back(DOTBRACKET__PAIRING_LEFT);					// '('
		}
		else if(j < i)
		{
			output.push_back(DOTBRACKET__PAIRING_RIGHT);				// ')'
		}
	}
	
	
	// Check if the number of parenthesis; '(' and ')' chars, is identical
#if DEBUG
	unsigned int left = 0;
	unsigned int right = 0;
	size_t k;
	
	for(k = 0; k < (size_t) n; k++)
	{
		switch(output[k])
		{
			case DOTBRACKET__PAIRING_LEFT:
				left++;
				break;
			case DOTBRACKET__PAIRING_RIGHT:
				right++;
				break;
			default:
				break;
		}
	}
	
	if(left != right)
	{
		throw std::out_of_range("DotBracket::format(): produced invalid 2D-structure '" + output + "', probably because of unexpected internal bug elsewhere in the code.");
	}
	
#endif // DEBUG
}



/**
 * @brief Matches a dotbracked string with a pattern in which unknown chars are encoded with a '?'
 *
 * @param dot_bracket_pattern A dotbracket string like "((..))((()))" which allows question tags for unknown structures: "(((???)))"
 * 
 * @todo make static?! it as no intric relation with the the classes content
 */
bool DotBracket::match(std::string &dot_bracket_pattern, std::string &dot_bracket_subject)
{
#if DEBUG
	if(dot_bracket_pattern.size() != dot_bracket_subject.size())
	{
		throw std::invalid_argument("DotBracket::match: comparing strings of different sizes (" + std::to_string((int) dot_bracket_pattern.size()) + ", " + std::to_string((int) dot_bracket_subject.size()) + ")");
	}
#endif //DEBUG
	
	for(std::string::size_type i = 0; i < dot_bracket_pattern.size(); ++i)
	{
		if(dot_bracket_pattern[i] != '?' && dot_bracket_pattern[i] != dot_bracket_subject[i])
		{
			return false;
		}
	}
	
	return true;
}
