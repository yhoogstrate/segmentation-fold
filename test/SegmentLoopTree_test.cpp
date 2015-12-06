/**
 * @file test/SegmentTree_test.cpp
 *
 * @date 2015-12-05
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */



#define BOOST_TEST_MODULE SegmentLoopTree


#include "main.hpp"

#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Position.hpp"
#include "SubSequence.hpp"
#include "Sequence.hpp"

#include "SegmentLoop.hpp"
#include "SegmentLoopTree.hpp"

#include <array>
#include <vector>

#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(Testing)



/**
 * @brief Tests quering using a SubSequence element instead of deep copied sub sequences
 *
 * @test
 *
 * @date 2015-12-05
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	std::string segmentloop_name = "Arteficial SegmentLoop construct";
	float energy = -1.234;
	
	Sequence sequence     = Sequence("ACUUGaGUA");
	std::vector <Pair> bonds = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop = SegmentLoop(segmentloop_name, sequence, bonds, energy);
	
	
	Sequence query_sequence     = Sequence("gACUUGaGUAg");
	
	Position p1a = query_sequence.data.begin() + 1;
	Position p1b = query_sequence.data.begin() + 10;
	
	BOOST_CHECK_EQUAL(*p1a , Nucleotide::A);
	BOOST_CHECK_EQUAL(*p1b , Nucleotide::G);
	
	SubSequence p1 = SubSequence(p1a, p1b);
	
	SegmentLoopTree segmentloops = SegmentLoopTree();
}

BOOST_AUTO_TEST_SUITE_END()
