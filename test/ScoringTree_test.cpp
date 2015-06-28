/**
 * @file test/ScoringTree_test.cpp
 *
 * @date 2014-04-20
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



#define BOOST_TEST_MODULE ScoringTree



#include "main.hpp"
#include "Utils/utils.hpp"

#include <array>

#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "Segment.hpp"
#include "ScoringTree.hpp"



#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests if a scoring tree with unsigned int's can insert and search items properly back
 *
 * @section DESCRIPTION
 * inserting and searching is in the same sequential, although the values are not initially sorted
 *
 * @test
 *
 * @date 2014-04-20
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	unsigned int value_1 =  5;
	unsigned int value_2 =  4;
	unsigned int value_3 =  1;
	unsigned int value_4 =  9;
	unsigned int value_5 = 67;
	
	Pair key_1 = Pair(5, 5);
	Pair key_2 = Pair(4, 4);
	Pair key_3 = Pair(1, 1);
	Pair key_4 = Pair(9, 9);
	Pair key_5 = Pair(6, 7);
	
	Pair key_6 = Pair(1, 9);
	
	ScoringTree<Pair, unsigned int>tree = ScoringTree<Pair, unsigned int>();
	tree.insert(key_1, value_1);
	tree.insert(key_2, value_2);
	tree.insert(key_3, value_3);
	tree.insert(key_4, value_4);
	tree.insert(key_5, value_5);
	
	unsigned int *value_back_1 = tree.search(key_1);
	unsigned int *value_back_2 = tree.search(key_2);
	unsigned int *value_back_3 = tree.search(key_3);
	unsigned int *value_back_4 = tree.search(key_4);
	unsigned int *value_back_5 = tree.search(key_5);
	unsigned int *value_back_6 = tree.search(key_6);
	
	BOOST_CHECK(*value_back_1 == value_1);
	BOOST_CHECK(*value_back_2 == value_2);
	BOOST_CHECK(*value_back_3 == value_3);
	BOOST_CHECK(*value_back_4 == value_4);
	BOOST_CHECK(*value_back_5 == value_5);
	BOOST_CHECK(value_back_6 == NULL);
}

/**
 * @brief Tests if overwriting is indeed disabled
 *
 * @test
 *
 * @date 2014-04-20
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	unsigned int value_1a = 5;
	unsigned int value_2a = 4;
	unsigned int value_1b = 50;
	unsigned int value_2b = 40;
	
	Pair key_1 = Pair(5, 5);
	Pair key_2 = Pair(4, 4);
	
	ScoringTree<Pair, unsigned int>tree = ScoringTree<Pair, unsigned int>();
	tree.insert(key_1, value_1a);
	tree.insert(key_2, value_2a);
	
#if DEBUG
	BOOST_CHECK_THROW(tree.insert(key_1, value_1b), std::invalid_argument);
	BOOST_CHECK_THROW(tree.insert(key_2, value_2b), std::invalid_argument);
#else //DEBUG
	tree.insert(key_1, value_1b);
	tree.insert(key_2, value_2b);
	
	unsigned int *value_back_1 = tree.search(key_1);
	unsigned int *value_back_2 = tree.search(key_2);
	
	BOOST_CHECK(*value_back_1 != value_1b);
	BOOST_CHECK(*value_back_2 != value_2b);
	BOOST_CHECK(*value_back_1 == value_1a);
	BOOST_CHECK(*value_back_2 == value_2a);
#endif //DEBUG
}

/**
 * @brief tests if inserting a segment based on a pair succeeds
 *
 * @test
 *
 * @date 2014-04-20
 *
 * @todo check whether the segment is a deep-copy - otherwise make the test with *Segment
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	Pair key_1 = Pair(5, 5);
	Pair key_2 = Pair(4, 4);
	
	ScoringTree<Pair, Segment>tree = ScoringTree<Pair, Segment>();
	
	std::string segment_name = "C/D-box K-turn";
	Sequence sequence_5p = Sequence("ACUUG");
	std::vector <Pair> bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence sequence_3p = Sequence("AUG");
	float energy = -1.234;
	Segment segment = Segment(segment_name, sequence_5p, bonds, sequence_3p, energy);
	
	tree.insert(key_1, segment);
	
	Segment *segment_1 = tree.search(key_1);
	Segment *segment_2 = tree.search(key_2);
	
	BOOST_CHECK(segment_1 != NULL && segment_1->get_gibbs_free_energy() == energy);
	BOOST_CHECK(segment_2 == NULL);
}
BOOST_AUTO_TEST_SUITE_END()
