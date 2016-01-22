/**
 * @file test/SegmentTraceback_test.cpp
 *
 * @date 2016-01-21
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
 * @date 2016-01-21
 */
BOOST_AUTO_TEST_CASE(Test_1)
{
	std::vector <Pair> bonds = {{Pair({0, 0}), Pair({0, 0}), Pair({0, 0})}};
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 3);
	
	unsigned int i = 0;
	unsigned int j = 20;
	
	while(segmenttraceback.traceback(i, j))
	{
		BOOST_CHECK_EQUAL(segmenttraceback.size() , 3);
	}
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 3);
}



/**
 * @brief Tests Segment::size() and whether it doesn't get affected by popping
 *
 * @test Segment::size()
 *
 * @date 2016-01-21
 */
BOOST_AUTO_TEST_CASE(Test_2)
{
	std::vector <Pair> bonds = {};
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 0);
	
	unsigned int i = 0;
	unsigned int j = 0;
	
	while(segmenttraceback.traceback(i, j))
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
 * @date 2016-01-21
 */
BOOST_AUTO_TEST_CASE(Test_3)
{
	std::vector <Pair> bonds = {{Pair({1, 1}), Pair({1, 1}), Pair({1, 1}), Pair({1, 1}), Pair({1, 1}), Pair({1, 1}), Pair({1, 1}), Pair({1, 1}), Pair({1, 1})}};
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 9);
	
	unsigned int i = 0;
	unsigned int j = 8;
	
	while(segmenttraceback.traceback(i, j))
	{
		BOOST_CHECK_EQUAL(segmenttraceback.size() , 9);
	}
	
	BOOST_CHECK_EQUAL(segmenttraceback.size() , 9);
}



/**
 * @brief Tests the Segment::pop() function
 *
 * @date 2016-01-21
 *
 * @test Segment::traceback()
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
 */
BOOST_AUTO_TEST_CASE(Test_4)
{
	// [1,1] , [2,3] -> differentiate -> [1,1] -> [1,2]
	std::vector <Pair> bonds = {{Pair({1, 1}), Pair({1, 2}) }};
	
	unsigned int i = 12;
	unsigned int j = 27;
	
	SegmentTraceback segmenttraceback = SegmentTraceback(bonds);
	
	BOOST_CHECK(segmenttraceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 13);//12 + 1 = 13
	BOOST_CHECK_EQUAL(j, 26);//27 - 1 = 26
	
	BOOST_CHECK(segmenttraceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 14);//12 + 1 + 1 = 14
	BOOST_CHECK_EQUAL(j, 24);//27 - 1 - 2 = 24
	
	BOOST_CHECK(segmenttraceback.traceback(i, j) == false);
	BOOST_CHECK_EQUAL(i, 14);//12 + 1 + 1 + 0 = 15
	BOOST_CHECK_EQUAL(j, 24);//23 - 1 - 2 - 0 = 24
	
	BOOST_CHECK(segmenttraceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 15);//12 + 1 + 1 + 0 + 1 = 15
	BOOST_CHECK_EQUAL(j, 23);//23 - 1 - 2 - 0 - 1 = 23
	
	BOOST_CHECK(segmenttraceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 16);//12 + 1 + 1 + 0 + 1 + 1 = 16
	BOOST_CHECK_EQUAL(j, 21);//23 - 1 - 2 - 0 - 1 - 2 = 21
	
	BOOST_CHECK(segmenttraceback.traceback(i, j) == false);
	BOOST_CHECK_EQUAL(i, 16);//12 + 1 + 1 + 0 + 1 + 1 + 0 = 16
	BOOST_CHECK_EQUAL(j, 21);//23 - 1 - 2 - 0 - 1 - 2 - 0 = 21
}

BOOST_AUTO_TEST_SUITE_END()
