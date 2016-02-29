/**
 * @file test/SegmentLoop_test.cpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2016 Youri Hoogstrate
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
 * @section DESCRIPTION
 * <PRE>
 *    Alignment:  Bonds relative to previous pair:  Traceback relative to (i,j):
 * i) ACUUG\      +1,-1                             +1,-1
 *    | | | a     +2,-1                             +3,-2
 * j) A U G/      +2,-1                             +5,-3
 * </PRE>
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	std::string segmentloop_name = "Arteficial SegmentLoop construct";
	float energy = -1.234f;
	
	Sequence sequence     = Sequence("ACUUGaGUA");
	std::vector <Pair> bonds = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop = SegmentLoop(segmentloop_name, sequence, bonds, energy);
	
	BOOST_CHECK_EQUAL(segmentloop.size(), sequence.size());
	BOOST_CHECK_EQUAL(segmentloop_name.compare(segmentloop.name),  0);
	
	unsigned int i = 0;
	unsigned int j = 1000;
	
	BOOST_CHECK(segmentloop.traceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 0    + 1);
	BOOST_CHECK_EQUAL(j, 1000 - 1);
	
	BOOST_CHECK(segmentloop.traceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 0    + 1 + 2);
	BOOST_CHECK_EQUAL(j, 1000 - 1 - 1);
	
	BOOST_CHECK(segmentloop.traceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 0    + 1 + 2 + 2);
	BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1);
	
	BOOST_CHECK(segmentloop.traceback.traceback(i, j) == false);
	BOOST_CHECK_EQUAL(i, 0    + 1 + 2 + 2 + 0);
	BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0);
	
	
	// Second (but identical) traceback; checks if reset function works properly
	BOOST_CHECK(segmentloop.traceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 0    + 1 + 2 + 2 + 0   + 1);
	BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1);
	
	BOOST_CHECK(segmentloop.traceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 0    + 1 + 2 + 2 + 0   + 1 + 2);
	BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1);
	
	BOOST_CHECK(segmentloop.traceback.traceback(i, j));
	BOOST_CHECK_EQUAL(i, 0    + 1 + 2 + 2 + 0   + 1 + 2 + 2);
	BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1 - 1);
	
	BOOST_CHECK(segmentloop.traceback.traceback(i, j) == false);
	BOOST_CHECK_EQUAL(i, 0    + 1 + 2 + 2 + 0   + 1 + 2 + 2 + 0);
	BOOST_CHECK_EQUAL(j, 1000 - 1 - 1 - 1 - 0   - 1 - 1 - 1 + 0);
}

BOOST_AUTO_TEST_SUITE_END()
