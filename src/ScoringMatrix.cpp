/**
 * @file src/ScoringMatrix.cpp
 *
 * @date 2016-01-20
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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
 * </PRE>
 */



#include "Pair.hpp"

#include "Nucleotide.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"
#include "ScoringMatrix.hpp"

#include "ZukerTraceback.hpp"



/**
 * @brief
 *
 * @date 2015-06-26
 *
 * @todo work with templates;
 * - Float for Energy values (V & W matrix)
 * - Int for identifiers (Segment objects)
 *
 * @section DESCRIPTION
 * For a matrix of n*n, all values (x,x) for 0 < x < n are initialized with a: 0
 * Similarly, all values (x,x-1) for 1 < x < n are initialized with a: 0
 *
 * Thus, only the positions in the matrix where x > y have to be located.
 * For a matrix of size: 4 * 4 (i = initation, s = store):
 *
 * <PRE>
 * [0,0][1,0][2,0][3,0]   [ i ][ s ][ s ][ s ]
 * [0,1][1,1][2,1][3,1]   [ i ][ i ][ s ][ s ]
 * [0,2][1,2][2,2][3,2]   [   ][ i ][ i ][ s ]
 * [0,3][1,3][2,3][3,3]   [   ][   ][ i ][ i ]
 * </PRE>
 *
 * The number of 's' type elements can be found using:
 * \f$s = \frac{n-1}{2} \times n\f$
 *
 * <PRE>
 * [ i ][ s ][ s ][ s ][ s ]
 * [ i ][ i ][ s ][ s ][ s ]
 * [   ][ i ][ i ][ s ][ s ]
 * [   ][   ][ i ][ i ][ s ]
 * [   ][   ][   ][ i ][ i ]
 * </PRE>
 */
template <class T>
ScoringMatrix<T>::ScoringMatrix(size_t arg_length, T arg_initialization_value):
	m(this->number_of_elements(arg_length))
{
	this->grid_size = arg_length;///@todo double check whether grid_size shouldn't be renamed to sequence_size or 1D size
	this->initialization_value = arg_initialization_value;
}



/**
 * @brief Gets the value at matrix point [x,y]
 *
 * @date 2016-01-20
 */
template <class T>
T ScoringMatrix<T>::get(Pair &pair)
{
	signed int position = this->get_position(pair);
	
	return (position < 0) ? this->initialization_value : this->m[(size_t) position];
}



/**
 * @brief Calculates the position (x,y) in the matrix corresponds to the position in the vector
 *
 * @date 2016-01-20
 *
 * @todo MAKE CONSTANT OF -2
 * @todo MAKE CONSTANT OF -1
 *
 * @todo put constant in ScoringMatrix.hpp
 */
template <class T>
signed int ScoringMatrix<T>::get_position(Pair &p)
{
#if DEBUG
	signed int output;
	
	if(p.first < this->grid_size && p.second < this->grid_size)
	{
		if(p.second == p.first || (p.first - 1) == p.second)			// The first 2 diagonals get -1 because they use a pre-defined variable
		{
			return -1;
		}
		else if(p.second > p.first)
		{
			size_t r = this->grid_size - p.first - 1;
			size_t bottom = this->number_of_elements(r);
			size_t row = (this->grid_size - 1) - p.second;
			
			///@todo change tis into size_t, and if out of bounds have to be returned, change this to MAX VALsize_t -1 and MAX VAL size_t -2
			return (signed int) this->m.size() - (signed int) bottom - (signed int) row - 1;
		}
		
	}
	
	// The out of bound may never occur in a good implementation of the algorithm
	throw std::invalid_argument("ScoringMatrix::get_position: Out of bound");
	
	return output;
	
#else //DEBUG
	
	// This piece also returns if a position in the upper triangle is requested, while it should throw out of bound
	// This exception will only be thrown if there is compiled with DEBUG
	
	///@todo reduce to one comparison instead of 2
	///@todo create function is_within_first_diagonals()
	if(p.second == p.first || (p.first - 1) == p.second)
	{
		return -1;
	}
	else
	{
		unsigned int r = this->grid_size - p.first - 1;
		unsigned int bottom = this->number_of_elements(r);
		unsigned int row = (this->grid_size - 1) - p.second;
	
		return this->m.size() - bottom - row - 1;
	}
# endif // DEBUG
}



/**
 * @brief Sets a value in the matrix
 *
 * @date 2015-12-09
 *
 * @todo Check whether diagonals are initiated; it takes unnecessairy computations
 */
template <class T>
void ScoringMatrix<T>::set(Pair &pair, T arg_value)
{
#if DEBUG
	signed int position = this->get_position(pair);
	
	if(position >= 0)
	{
		this->m[(size_t) position] = arg_value;
	}
	else
	{
		throw std::invalid_argument("ScoringMatrix::set: Out of bound");
	}
#else //DEBUG
	this->m[(size_t) this->get_position(pair)] = arg_value;
#endif //DEBUG
}



/**
 * @brief Returns the number of reserved elements in the vector
 *
 * @date 2015-06-24
 */
template <class T>
size_t ScoringMatrix<T>::size(void)
{
	return this->m.size();
}



/**
 * @brief Calculates the number of elements in the array that correspond to the size of the grid
 *
 * @section DESCRIPTION
 * Finds the number of 's' elements in a matrix of size n*n:
 * <PRE>
 * [ i ][ s ][ s ][ s ][ s ]
 * [ i ][ i ][ s ][ s ][ s ]
 * [   ][ i ][ i ][ s ][ s ]
 * [   ][   ][ i ][ i ][ s ]
 * [   ][   ][   ][ i ][ i ]
 * </PRE>
 *
 * Using the formula:
 * \f$s = \frac{n-1}{2} \times n\f$
 *
 * @date 2013-09-19
 *
 * @note "((n - 1) * n) / 2" works, while "((n - 1) / 2) * n" does not because of floating point divisions
 */
template <class T>
size_t ScoringMatrix<T>::number_of_elements(size_t n)
{
	return ((n - 1) * n) / 2;
}



/**
 * @brief Fill the remainder of the matrix with the init value.
 *
 * @todo use m.assign()
 */
template <class T>
void ScoringMatrix<T>::fill(T arg_value)
{
	unsigned int i;
	
	for(i = 0; i < this->m.size(); i++)
	{
		this->m[i] = arg_value;
	}
}


// Explicit instantiation:
template class ScoringMatrix<signed int>;
template class ScoringMatrix<float>;
template class ScoringMatrix<SegmentTraceback *>;
template class ScoringMatrix<char>;
template class ScoringMatrix<Pair>;

template class ScoringMatrix<traceback_jump2>;

//template class ScoringMatrix<unsigned int>;
//template class ScoringMatrix<double>;
//template class ScoringMatrix<short>;
//template class ScoringMatrix<bool>; <- this one causes problems! do never use it.
