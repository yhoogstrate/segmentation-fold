/**
 * @file test/SegmentTraceback_test.cpp
 *
 * @date 2015-05-02
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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
 * <PRE>
 */



#define BOOST_TEST_MODULE SegmentTraceback


#include "main.hpp"
#include <array>


#include "Pair.hpp"
#include "SegmentTraceback.hpp"



#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests the Segment::size() function
 *
 * @test Tests Segment::size() and whether it doesn't get affected by popping
 *
 * @date 2015-08-05
 */
BOOST_AUTO_TEST_CASE(Test_size_1)
{
	std::vector <Pair> bonds = {{Pair({0, 0}), Pair({2, 1}), Pair({4, 2})}};
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 3);
	
	signed int i;
	signed int j;
	
	while(segmenttraceback.pop(i, j))
	{
		BOOST_CHECK_EQUAL(segmenttraceback.size() , 3);
	}
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 3);
}



/**
 * @brief Tests Segment::size() and whether it doesn't get affected by popping
 *
 * @test Segment::size()
 *  *
 * @date 2015-08-05
 */
BOOST_AUTO_TEST_CASE(Test_size_2)
{
	std::vector <Pair> bonds = {};
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 0);
	
	signed int i;
	signed int j;
	
	while(segmenttraceback.pop(i, j))
	{
		BOOST_REQUIRE_MESSAGE(1 != 1, "An empty segmenttraceback should not pop");
	}
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 0);
}



/**
 * @brief Tests Segment::size() and whether it doesn't get affected by popping
 *
 * @test Segment::size()
 *
 * @date 2015-08-05
 */
BOOST_AUTO_TEST_CASE(Test_size_3)
{
	std::vector <Pair> bonds = {{Pair({0, 0}), Pair({1, 1}), Pair({2, 2}), Pair({3, 3}), Pair({4, 4}), Pair({5, 5}), Pair({6, 6}), Pair({7, 7}), Pair({8, 8})}};
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 9);
	
	signed int i;
	signed int j;
	
	while(segmenttraceback.pop(i, j))
	{
		BOOST_CHECK_EQUAL(segmenttraceback.size() , 9);
	}
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 9);
}



/**
 * @brief Tests the Segment::pop() function
 * 
 * @test Segment::pop()
 * 
 * @section DESCRIPTION
 * The following structure is an example:
 * 
 * <PRE>
 * 5') ............(i). ..(i')..
 *       |||  |||   | : :  |    .
 * 3') ............(j)....(j')..
 * </PRE>
 *
 * If the bonds are parsed correctly, it would be a vector like this:
 * 
 * [ [1, 1], [2, 3] ], because the first bond after i,j is [i+1, j-1] and the second bond is [i+2, j-3]
 * 
 * @date 2015-08-05
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	std::vector <Pair> bonds = {{Pair({1, 1}), Pair({2, 3}) }};
	
	int i = 12;
	int j = 27;
	
	int i_add, j_sub, i_add_tmp, j_sub_tmp;
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	
	// Second (but identical) traceback; checks if reset function works properly
	BOOST_CHECK(segmenttraceback.pop(i_add, j_sub));
	BOOST_CHECK_EQUAL(i_add, 1);
	BOOST_CHECK_EQUAL(j_sub, 1);
	BOOST_CHECK_EQUAL(i+i_add, 13);
	BOOST_CHECK_EQUAL(j-j_sub, 26);
	
	BOOST_CHECK(segmenttraceback.pop(i_add, j_sub));
	BOOST_CHECK_EQUAL(i_add, 2);
	BOOST_CHECK_EQUAL(j_sub, 3);
	BOOST_CHECK_EQUAL(i+i_add, 14);
	BOOST_CHECK_EQUAL(j-j_sub, 24);
	
	i_add_tmp = i_add;
	j_sub_tmp = j_sub;
	
	BOOST_CHECK(segmenttraceback.pop(i_add, j_sub) == false);
	BOOST_CHECK_EQUAL(i_add_tmp , i_add);// Check that popping a false doesn't change the values
	BOOST_CHECK_EQUAL(j_sub_tmp , j_sub);
}



/**
 * @brief Tests the Segment::pop() function
 * 
 * @test Segment::pop_traceback()
 * 
 * @section DESCRIPTION
 * The following structure is an example:
 * 
 * <PRE>
 * 5') ............(i). ..(i')..
 *       |||  |||   | : :  |    .
 * 3') ............(j)....(j')..
 * </PRE>
 *
 * If the bonds are parsed correctly, it would be a vector like this:
 * 
 * [ [1, 1], [2, 3] ], because the first bond after i,j is [i+1, j-1] and the second bond is [i+2, j-3]
 * 
 * @date 2015-08-05
 */
BOOST_AUTO_TEST_CASE(Test2b)
{
	std::vector <Pair> bonds = {{Pair({1, 1}), Pair({2, 3}) }};
	
	unsigned int i = 12;
	unsigned int j = 27;
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK(segmenttraceback.pop_traceback(i, j));
	BOOST_CHECK_EQUAL(i, 13);//12 + 1 = 13
	BOOST_CHECK_EQUAL(j, 26);//27 - 1 = 26
	
	BOOST_CHECK(segmenttraceback.pop_traceback(i, j));
	BOOST_CHECK_EQUAL(i, 15);//13 + 2 = 15
	BOOST_CHECK_EQUAL(j, 23);//26 - 3 = 23
	
	BOOST_CHECK(segmenttraceback.pop_traceback(i, j) == false);
	BOOST_CHECK_EQUAL(i, 15);//15 + 0 = 15
	BOOST_CHECK_EQUAL(j, 23);//23 - 0 = 23
	
	BOOST_CHECK(segmenttraceback.pop_traceback(i, j));
	BOOST_CHECK_EQUAL(i, 16);//15 + 1 = 16
	BOOST_CHECK_EQUAL(j, 22);//23 - 1 = 22
	
	BOOST_CHECK(segmenttraceback.pop_traceback(i, j));
	BOOST_CHECK_EQUAL(i, 18);//16 - 2
	BOOST_CHECK_EQUAL(j, 19);//22 - 3 = 19
	
	BOOST_CHECK(segmenttraceback.pop_traceback(i, j) == false);
	BOOST_CHECK_EQUAL(i, 18);
	BOOST_CHECK_EQUAL(j, 19);//19 - 0
}



/**
 * @brief Tests the Segment::pop() and Segment::reset() functions
 *
 * @test Segment::pop Segment::reset()
 *
 * @date 2015-04-22
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	int i, j, i_tmp, j_tmp;
	
	std::vector <Pair> bonds = {{Pair({1, 2}), Pair({2, 4}), Pair({3, 6})}};
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	unsigned int k;
	for(k = 1; k < 25; k++)
	{
		BOOST_CHECK(segmenttraceback.pop(i, j));
		BOOST_CHECK_EQUAL(i , 1);
		BOOST_CHECK_EQUAL(j , 2);
		
		BOOST_CHECK(segmenttraceback.pop(i, j));
		BOOST_CHECK_EQUAL(i , 2);
		BOOST_CHECK_EQUAL(j , 4);
		
		BOOST_CHECK(segmenttraceback.pop(i, j));
		BOOST_CHECK_EQUAL(i , 3);
		BOOST_CHECK_EQUAL(j , 6);
		
		i_tmp = i;
		j_tmp = j;
		
		BOOST_CHECK(segmenttraceback.pop(i, j) == false);
		BOOST_CHECK_EQUAL(i_tmp , i);// Check that popping a false doesn't change the values:
		BOOST_CHECK_EQUAL(j_tmp , j);
	}
}

BOOST_AUTO_TEST_SUITE_END()
