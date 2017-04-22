/**
 * @file include/ScoringMatrix.hpp
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



#ifndef SCORINGMATRIX_HPP
#define	SCORINGMATRIX_HPP



#include "main.hpp"


/**
 * @brief A (memory efficient) lower triangle of an n*n scoring matrix plus an additional vertical.
 *
 * @section DESCRIPTION
 * For a matrix of n*n, all values (x,x) for 0 < x < n are initialized with a: 0
 * Similarly, all values (x,x-1) for 1 < x < n are initialized with a: 0
 *
 * Thus, only the positions in the matrix where x > y have to be located.
 * For a matrix of size: 4 * 4 (i = NOT_YET_CALCULATED, s = un initialized):

[0,0][1,0][2,0][3,0]   [ i ][ s ][ s ][ s ]
[0,1][1,1][2,1][3,1]   [ i ][ i ][ s ][ s ]
[0,2][1,2][2,2][3,2]   [   ][ i ][ i ][ s ]
[0,3][1,3][2,3][3,3]   [   ][   ][ i ][ i ]

//size = ((n-1)/2) * n + n + (n-1) = ((((n-1)/2)+2) * n) - 1 = (((n+3)/2) * n) -1
//size = max(size,0)

[ i ][ s ][ s ][ s ][ s ]
[ i ][ i ][ s ][ s ][ s ]
[   ][ i ][ i ][ s ][ s ]
[   ][   ][ i ][ i ][ s ]
[   ][   ][   ][ i ][ i ]
 */
template <class T>
class ScoringMatrix
{
	private:
		size_t grid_size;
		
		T initialization_value;
		std::vector<T> m;
		
	public:
		ScoringMatrix(size_t arg_length, T arg_initialization_value);///@todo use size_t instead of unsigned int?
		
		signed int get_position(Pair &p);
		
		T get(Pair &p);
		
		void set(Pair &p, T arg_value);
		
		inline size_t size(void);
		static size_t number_of_elements(size_t n);// static, as it uses n as argument instead of this->size()
		
		///@todo Nucleotide operator[](size_t); << directly obtain from the vector this->m, allows caching of its iterator
		
		void fill(T arg_value);
};

#endif	// SCORINGMATRIX_HPP
