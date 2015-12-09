/**
 * @file test/ScoringMatrix_test.cpp
 *
 * @date 2015-12-09
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
 * along with this program. if not, see <http://www.gnu.org/licenses/>.
 */



#define BOOST_TEST_MODULE ScoringMatrix

#include <boost/test/included/unit_test.hpp>


#include "Pair.hpp"

#include "Nucleotide.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"

#include "ScoringMatrix.hpp"



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether the size of the matrix scales correctly up with the size of a Sequence
 *
 * @test
 *
 * @date 2015-06-24
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	unsigned int n_elements;
	
	for(unsigned int i = 0; i <= 120; i++)
	{
		ScoringMatrix<signed int> matrix = ScoringMatrix<signed int>(i, 0);
		
		switch(i)
		{
			case 1:
				n_elements = 0;
				
				BOOST_CHECK_EQUAL(matrix.size(), n_elements);
				BOOST_CHECK_EQUAL(matrix.number_of_elements(i), n_elements);
				break;
			case 2:
				n_elements = 1;
				
				BOOST_CHECK_EQUAL(matrix.size(), n_elements);
				BOOST_CHECK_EQUAL(matrix.number_of_elements(i), n_elements);
				break;
			case 3:
				n_elements = 3;
				
				BOOST_CHECK_EQUAL(matrix.size(), n_elements);
				BOOST_CHECK_EQUAL(matrix.number_of_elements(i), n_elements);
				
				BOOST_CHECK_EQUAL(matrix.size(), n_elements);
				BOOST_CHECK_EQUAL(matrix.number_of_elements(i), n_elements);
				break;
			case 4:
				n_elements = 6;
				
				BOOST_CHECK_EQUAL(matrix.size(), n_elements);
				BOOST_CHECK_EQUAL(matrix.number_of_elements(i), n_elements);
				break;
			case 5:
				n_elements = 10;
				
				BOOST_CHECK_EQUAL(matrix.size(), n_elements);
				BOOST_CHECK_EQUAL(matrix.number_of_elements(i), n_elements);
				break;
			case 6:
				n_elements = 15;
				
				BOOST_CHECK_EQUAL(matrix.size(), n_elements);
				BOOST_CHECK_EQUAL(matrix.number_of_elements(i), n_elements);
				break;
			default:
				BOOST_CHECK_EQUAL(matrix.size(), ((i - 1) * i) / 2);
				BOOST_CHECK_EQUAL(matrix.size(), matrix.number_of_elements(i));
				break;
		}
	}
}



/**
 * @brief Tests whether the correct element id's of the vector are returned
 *
 * @test
 *
 * @date 2015-12-09
 *
 * @note This test is testing a function that originally was planned to be private. Hopefully one day we succeed in testing a private function: http://boost.2283326.n4.nabble.com/Testing-private-methods-in-Boost-td2599819.html
 *
 * Positions:
 *
 * <PRE>
	[-1 ][-1 ][-2 ][-2 ][-2 ]
	[ 0 ][-1 ][-1 ][-2 ][-2 ]
	[ 1 ][ 4 ][-1 ][-1 ][-2 ]
	[ 2 ][ 5 ][ 7 ][-1 ][-1 ]
	[ 3 ][ 6 ][ 8 ][ 9 ][-1 ]
 * </PRE>
 *
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	ScoringMatrix<signed int> matrix = ScoringMatrix<signed int>(5, 0);
	Pair p = Pair();
	
	// Test correct elements from the lower triangle:
	p.first = 0;
	p.second = 1;
	BOOST_CHECK_EQUAL(0, matrix.get_position(p));
	p.second = 2;
	BOOST_CHECK_EQUAL(1, matrix.get_position(p));
	p.second = 3;
	BOOST_CHECK_EQUAL(2, matrix.get_position(p));
	p.second = 4;
	BOOST_CHECK_EQUAL(3, matrix.get_position(p));
	
	p.first = 1;
	p.second = 2;
	BOOST_CHECK_EQUAL(4, matrix.get_position(p));
	p.second = 3;
	BOOST_CHECK_EQUAL(5, matrix.get_position(p));
	p.second = 4;
	BOOST_CHECK_EQUAL(6, matrix.get_position(p));
	
	p.first = 2;
	p.second = 3;
	BOOST_CHECK_EQUAL(7, matrix.get_position(p));
	p.second = 4;
	BOOST_CHECK_EQUAL(8, matrix.get_position(p));
	
	p.first = 3;
	p.second = 4;
	BOOST_CHECK_EQUAL(9, matrix.get_position(p));
	
	
	// Test correct element positions of the diagonals
	for(p.first = 0; p.first <= 4; p.first++)
	{
		p.second = p.first;
		BOOST_CHECK_EQUAL(-1, matrix.get_position(p));
	}
	
	// Test correct element positions of the diagonals
	for(p.first = 1; p.first <= 4; p.first++)
	{
		p.second = p.first - 1;
		BOOST_CHECK_EQUAL(-1, matrix.get_position(p));
	}
	
	
#if DEBUG
	// Located at the not allocated upper triangle:
	p.first = 2;
	p.second = 0;
	BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
	p.first = 3;
	p.second = 1;
	BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
	p.first = 4;
	p.second = 2;
	BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
	
	p.first = 3;
	p.second = 0;
	BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
	p.first = 4;
	p.second = 1;
	BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
	
	p.first = 4;
	p.second = 0;
	BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
	
	// Positions that fall outside the entire matrix:
	for(int i = 5; i <= 100; i += 5)
	{
		p.first = 2;
		p.second = i;
		BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
		p.first = i;
		p.second = 2;
		BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
		p.first = 5;
		p.second = i;
		BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
		p.first = i;
		p.second = 5;
		BOOST_CHECK_THROW(matrix.get_position(p), std::invalid_argument);
	}
#endif //DEBUG
}



/**
 * @brief Tests whether setting and getting works correctly
 *
 * @test
 *
 * @date 2015-12-09
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	unsigned int n, k;
	
	k = 0;
	
	Pair pair = Pair(0, 0);
	
	for(n = 0; n <= 150; n += 10)
	{
		ScoringMatrix<signed int> matrix = ScoringMatrix<signed int>(n, -1);
		
		for(pair.second = 0; pair.second < n; pair.second++)
		{
			for(pair.first = 0; pair.first <= pair.second; pair.first++)
			{
				if(pair.second == pair.first)							// Diagonal line, must return default value
				{
					BOOST_CHECK_EQUAL(matrix.get(pair), -1);
				}
				else													// Matrix cell, must return the requested value (k)
				{
					matrix.set(pair, k);
					BOOST_CHECK_EQUAL(matrix.get(pair), k++);
				}
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
