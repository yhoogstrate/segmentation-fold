/**
 * @file test/SegmentLoop_test.cpp
 *
 * @date 2015-12-03
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
 * </PRE>
 */



#define BOOST_TEST_MODULE SegmentLoop



#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "SegmentLoop.hpp"

#include <array>



#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests the SegmentLoop->size() function
 *
 * @test
 *
 * @date 2014-03-29
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	std::string segment_name = "Arteficial segmentloop construct";
	
	// ACUUG\
	// | | | a
	// A U G/
	
	/*
	Sequence sequence_5p = Sequence("ACUUGaGUC");
	std::vector <Pair> bonds = {{Pair({0, 2}), Pair({2, 1}), Pair({4, 0})}};///@todo incorrect bonds
	Sequence sequence_3p = Sequence("AUG");
	*/
	float energy = -1.234;
	/*
	SegmentLoop segment = SegmentLoop(segment_name, sequence_5p, bonds, sequence_3p, energy);
	
	BOOST_CHECK_EQUAL(segment.size(Direction::FivePrime),  5);
	BOOST_CHECK_EQUAL(segment.size(Direction::ThreePrime), 3);
	
	BOOST_CHECK_EQUAL(segment_name.compare(segment.name),  0);
	*/
}

///**
//* @brief Tests the SegmentLoop->pop() function
//*
//* @test
//*
//* @date 2015-08-06
//*/
//BOOST_AUTO_TEST_CASE(Test2)
//{
///*
//'Simple' case; linear traceback

//Alignment:      Bonds:      Traceback:
//UGUGAU          3-2         -1,-1
//|||          4-1         -2,-2
//AGU          5-0         -3,-3
//*/

//std::string segment_name   = "C/D-box K-turn";

//Sequence sequence_5p     = Sequence("UGUGAU");
//std::vector <Pair> bonds = {{Pair({4, 1}), Pair({1, 1}), Pair({1, 1})}};
//Sequence     sequence_3p = Sequence("UGA");

//float energy = -1.234;

//int i = 0;
//int j = 1000;

//SegmentLoop segment = SegmentLoop(segment_name, sequence_5p, bonds, sequence_3p, energy);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0 + 4);
//BOOST_CHECK_EQUAL(j, 1000 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 4 + 1);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 4 + 1 + 1);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j) == false);
//BOOST_CHECK_EQUAL(i, 0    + 4 + 1 + 1 + 0);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0);


//// Second (but identical) traceback; checks if reset function works properly
//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 4 + 1 + 1 + 0   + 4);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 4 + 1 + 1 + 0   + 4 + 1);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 4 + 1 + 1 + 0   + 4 + 1 + 1);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j) == false);
//BOOST_CHECK_EQUAL(i, 0    + 4 + 1 + 1 + 0   + 4 + 1 + 1 + 0);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1 - 1 + 0);
//}



///**
//* @brief Tests the SegmentLoop->pop() function
//*
//* @test
//*
//* @date 2015-08-06
//*/
//BOOST_AUTO_TEST_CASE(Test3)
//{
///*
//Alignment:      Bonds:      Traceback:
//AAAAAA          2-2         3,1
//| ||          4-1         2,1
//A AA          5-0         1,1
//*/

//std::string segment_name = "Articificial example";

//Sequence sequence_5p = Sequence("AAAAAA");
//std::vector <Pair> bonds = {{Pair({3, 1}), Pair({2, 1}), Pair({1, 1})}};
//Sequence sequence_3p = Sequence("AAA");

//float energy = -1.234;

//int i = 0;
//int j = 1000;

//SegmentLoop segment = SegmentLoop(segment_name, sequence_5p, bonds, sequence_3p, energy);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 3);
//BOOST_CHECK_EQUAL(j, 1000 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 3 + 2);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 3 + 2 + 1);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j) == false);
//BOOST_CHECK_EQUAL(i, 0    + 3 + 2 + 1 + 0);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0);


//// Second (but identical) traceback; checks if reset function works properly
//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 3 + 2 + 1 + 0   + 3);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 3 + 2 + 1 + 0   + 3 + 2);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j));
//BOOST_CHECK_EQUAL(i, 0    + 3 + 2 + 1 + 0   + 3 + 2 + 1);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1 - 1);

//BOOST_CHECK(segment.traceback.traceback(i, j) == false);
//BOOST_CHECK_EQUAL(i, 0    + 3 + 2 + 1 + 0   + 3 + 2 + 1 + 0);
//BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1 - 1 + 0);
//}


BOOST_AUTO_TEST_SUITE_END()
