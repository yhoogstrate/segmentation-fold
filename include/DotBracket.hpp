/**
 * @file include/DotBracket.hpp
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


#ifndef DOTBRACKET_HPP
#define	DOTBRACKET_HPP


#include "main.hpp"



/**
 * @brief Container of secondary structure as DotBracket format.
 *
 * @section DESCRIPTION
 * Maintains a DotBracket formatted 2D structure (without pseudo-knotting)
 *
 * Specification of the file format:
 * <a href="https://wiki.galaxyproject.org/Learn/Datatypes#Dbn">https://wiki.galaxyproject.org/Learn/Datatypes#Dbn</a>
 *
 * @date 2016-01-20
 */
class DotBracket
{
	private:
		std::vector <std::pair<signed int, signed int>> pairings;
		
	public:
		DotBracket();
		
		void store(unsigned int arg_i, unsigned int arg_j);
		int  find(unsigned int arg_i);
		void format(unsigned int n, std::string &output);
		
		bool match(std::string &dot_bracket_pattern, std::string &dot_bracket_subject);
};


#endif	// DOTBRACKET_HPP
