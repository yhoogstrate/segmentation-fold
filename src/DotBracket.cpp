/**
 * @file src/DotBracket.cpp
 *
 * @date 2016-01-20
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2015 Youri Hoogstrate
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
 * @date 2013-09-11
 *
 * @todo Make this a WHILE loop and avoid the early returns
 * @todo get rid of signed here - maybe make an at() function that requires to find something
*/
signed int DotBracket::find(unsigned int arg_i)
{
	unsigned int i;
	
	for(i = 0; i < this->pairings.size(); i++)
	{
		if(this->pairings[i].first == arg_i)
		{
			return (signed int) this->pairings[i].second;
		}
		else if(this->pairings[i].second == arg_i)
		{
			return (signed int) this->pairings[i].first;
		}
	}
	
	return UNBOUND;
}



/**
 * @brief Stores pairing (i,j)
 *
 * @param arg_i Nucleotide position in the RNA sequence pairing with (,j)
 * @param arg_j Nucleotide position in the RNA sequence pairing with (i,)
 *
 * @date 2013-09-11
*/
void DotBracket::store(unsigned int arg_i, unsigned int arg_j)
{
	this->pairings.push_back({arg_i, arg_j});
}



/**
 * @brief Formats and prints the 2D structure in DotBracket format
 *
 * @param rnaSequenceLength The length of the sequence of which the 2D strucutre is calculated
 *
 * @date 2016-01-20
 *
 * @todo Make it work with streams, so that output can be put in a file.
 * @todo Make the variables private const static class variables
 */
void DotBracket::format(unsigned int n, std::string &output)
{
	unsigned int i;
	signed int j;
	
	for(i = 0; i < n; i++)
	{
		j = this->find(i);
		
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
}



/**
 * @brief Matches a dotbracked string with a pattern in which unknown chars are encoded with a '?'
 *
 * @date 2015-04-20
 *
 * @param dot_bracket_pattern A dotbracket string like "((..))((()))" which allows question tags for unknown structures: "(((???)))"
 */
bool DotBracket::match(std::string &dot_bracket_pattern, std::string &dot_bracket_subject)
{
#if DEBUG

	if(dot_bracket_pattern.size() == dot_bracket_subject.size())
	{
		for(std::string::size_type i = 0; i < dot_bracket_pattern.size(); ++i)
		{
			if(dot_bracket_pattern[i] != '?' && dot_bracket_pattern[i] != dot_bracket_subject[i])
			{
				return false;
			}
		}
		
		return true;
	}
	else
	{
		return false;
	}
	
#else //DEBUG
	
	for(std::string::size_type i = 0; i < dot_bracket_pattern.size(); ++i)
	{
		if(dot_bracket_pattern[i] != '?' && dot_bracket_pattern[i] != dot_bracket_subject[i])
		{
			return false;
		}
	}
	
	return true;
	
#endif //DEBUG
}
