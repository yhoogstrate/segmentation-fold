/**
 * @file src/Nucleotide.hpp
 *
 * @date 20-apr-2014
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



/**
 * @brief Describes a DNA or RNA nucleotide
 *
 * @date 10-mar-2014
 *
 * @author Youri Hoogstrate
 *
 * @section DESCRIPTION
 * The type is casted to char to reduce the size of the memory, since
 * this is not automatically done by the compiler.
 */
enum struct Nucleotide : char
{
	A = 0,
	C = 1,
	G = 2,
	U, T = 3
};
