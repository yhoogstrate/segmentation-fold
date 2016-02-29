/**
 * @file test/SegmentTree_test.cpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
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
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	std::string segmentloop_name = "Arteficial SegmentLoop construct";
	float energy = -1.234f;
	
	Sequence sequence     = Sequence("ACUUGaGUA");
	std::vector <Pair> bonds = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop = SegmentLoop(segmentloop_name, sequence, bonds, energy);
	
	
	Sequence query_sequence     = Sequence("gACUUGaGUAg");
	
	Position p1a = query_sequence.data.begin() + 1;
	Position p1b = query_sequence.data.begin() + 1 + 8;
	
	BOOST_CHECK_EQUAL(*p1a , Nucleotide::A);
	BOOST_CHECK_EQUAL(*p1b , Nucleotide::A);
	
	SubSequence p1 = SubSequence(p1a, p1b);
	
	SegmentLoopTree segmentloop_tree = SegmentLoopTree();
	segmentloop_tree.insert(segmentloop);
	
	BOOST_CHECK(segmentloop_tree.search(p1) != nullptr);
	BOOST_CHECK(segmentloop_tree.search(p1)->sequence == sequence);
}



/**
 * @brief Tests whether left and right are chosen correctly
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	std::string segmentloop_name_01 = "A";
	float energy_01 = -1.00;
	Sequence sequence_01     = Sequence("AAAA");
	std::vector <Pair> bonds_01 = {{Pair({1, 1})}};
	SegmentLoop segmentloop_01 = SegmentLoop(segmentloop_name_01, sequence_01, bonds_01, energy_01);
	
	std::string segmentloop_name_02 = "B";
	float energy_02 = -2.00;
	Sequence sequence_02     = Sequence("AA");
	std::vector <Pair> bonds_02 = {{Pair({1, 1})}};
	SegmentLoop segmentloop_02 = SegmentLoop(segmentloop_name_02, sequence_02, bonds_02, energy_02);
	
	std::string segmentloop_name_03 = "C";
	float energy_03 = -3.00;
	Sequence sequence_03     = Sequence("AAA");
	std::vector <Pair> bonds_03 = {{Pair({1, 1})}};
	SegmentLoop segmentloop_03 = SegmentLoop(segmentloop_name_03, sequence_03, bonds_03, energy_03);
	
	
	//Requires the Test_SegmentLoopTree object because the SegmentLoopTree::root is private
	Test_SegmentLoopTree segmentloop_tree = Test_SegmentLoopTree();
	segmentloop_tree.insert(segmentloop_01);
	segmentloop_tree.insert(segmentloop_02);
	segmentloop_tree.insert(segmentloop_03);
	
	BOOST_CHECK(&(segmentloop_tree.root->segmentloop) == &segmentloop_01);
	
	BOOST_CHECK(&(segmentloop_tree.root->left->segmentloop) == &segmentloop_02);
	BOOST_CHECK(segmentloop_tree.root->right == nullptr);
	
	BOOST_CHECK(&(segmentloop_tree.root->left->right->segmentloop) == &segmentloop_03);
	BOOST_CHECK(segmentloop_tree.root->left->left == nullptr);
}



/**
 * @brief Tests whether the tree queries correctly when multiple elements are inserted
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	std::string segmentloop_name_01 = "A";
	float energy_01 = -1.00;
	Sequence sequence_01     = Sequence("AAAAA");
	std::vector <Pair> bonds_01 = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop_01 = SegmentLoop(segmentloop_name_01, sequence_01, bonds_01, energy_01);
	
	std::string segmentloop_name_02 = "B";
	float energy_02 = -2.00;
	Sequence sequence_02     = Sequence("CAAAA");
	std::vector <Pair> bonds_02 = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop_02 = SegmentLoop(segmentloop_name_02, sequence_02, bonds_02, energy_02);
	
	std::string segmentloop_name_03 = "C";
	float energy_03 = -3.00;
	Sequence sequence_03     = Sequence("ACAAA");
	std::vector <Pair> bonds_03 = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop_03 = SegmentLoop(segmentloop_name_03, sequence_03, bonds_03, energy_03);
	
	std::string segmentloop_name_04 = "D";
	float energy_04 = -4.00;
	Sequence sequence_04     = Sequence("AACAA");
	std::vector <Pair> bonds_04 = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop_04 = SegmentLoop(segmentloop_name_04, sequence_04, bonds_04, energy_04);
	
	std::string segmentloop_name_05 = "E";
	float energy_05 = -5.00;
	Sequence sequence_05     = Sequence("AACA");
	std::vector <Pair> bonds_05 = {{Pair({1, 1}), Pair({2, 1}) }};
	SegmentLoop segmentloop_05 = SegmentLoop(segmentloop_name_05, sequence_05, bonds_05, energy_05);
	
	std::string segmentloop_name_06 = "F";
	float energy_06 = -6.00;
	Sequence sequence_06     = Sequence("AAAAC");
	std::vector <Pair> bonds_06 = {{Pair({1, 1}), Pair({2, 1}), Pair({2, 1})}};
	SegmentLoop segmentloop_06 = SegmentLoop(segmentloop_name_06, sequence_06, bonds_06, energy_06);
	
	
	Sequence query_sequence     = Sequence("gAAAAAg");
	
	Position p1a = query_sequence.data.begin() + 1;
	Position p1b = query_sequence.data.begin() + 1 + 4;
	
	BOOST_CHECK_EQUAL(*p1a , Nucleotide::A);
	BOOST_CHECK_EQUAL(*p1b , Nucleotide::A);
	
	SubSequence p1 = SubSequence(p1a, p1b);
	
	BOOST_CHECK_EQUAL(p1.size , 5);
	
	SegmentLoopTree segmentloop_tree_01 = SegmentLoopTree();
	segmentloop_tree_01.insert(segmentloop_01);
	segmentloop_tree_01.insert(segmentloop_02);
	segmentloop_tree_01.insert(segmentloop_03);
	segmentloop_tree_01.insert(segmentloop_04);
	segmentloop_tree_01.insert(segmentloop_05);
	segmentloop_tree_01.insert(segmentloop_06);
	
	BOOST_CHECK(segmentloop_tree_01.search(p1) != nullptr);
	BOOST_CHECK(segmentloop_tree_01.search(p1)->gibbs_free_energy == -1.0);
	
	
	SegmentLoopTree segmentloop_tree_02 = SegmentLoopTree();
	segmentloop_tree_02.insert(segmentloop_06);
	segmentloop_tree_02.insert(segmentloop_05);
	segmentloop_tree_02.insert(segmentloop_04);
	segmentloop_tree_02.insert(segmentloop_03);
	segmentloop_tree_02.insert(segmentloop_02);
	segmentloop_tree_02.insert(segmentloop_01);
	
	BOOST_CHECK(segmentloop_tree_02.search(p1) != nullptr);
	BOOST_CHECK(segmentloop_tree_02.search(p1)->gibbs_free_energy == -1.0);
	
	
	SegmentLoopTree segmentloop_tree_03 = SegmentLoopTree();
	segmentloop_tree_03.insert(segmentloop_04);
	segmentloop_tree_03.insert(segmentloop_05);
	segmentloop_tree_03.insert(segmentloop_06);
	segmentloop_tree_03.insert(segmentloop_01);
	segmentloop_tree_03.insert(segmentloop_02);
	segmentloop_tree_03.insert(segmentloop_03);
	
	BOOST_CHECK(segmentloop_tree_03.search(p1) != nullptr);
	BOOST_CHECK(segmentloop_tree_03.search(p1)->gibbs_free_energy == -1.0);
	
	
	SegmentLoopTree segmentloop_tree_04 = SegmentLoopTree();
	segmentloop_tree_04.insert(segmentloop_01);
	segmentloop_tree_04.insert(segmentloop_06);
	segmentloop_tree_04.insert(segmentloop_02);
	segmentloop_tree_04.insert(segmentloop_05);
	segmentloop_tree_04.insert(segmentloop_03);
	segmentloop_tree_04.insert(segmentloop_04);
	
	BOOST_CHECK(segmentloop_tree_04.search(p1) != nullptr);
	BOOST_CHECK(segmentloop_tree_04.search(p1)->gibbs_free_energy == -1.0);
	
	
	SegmentLoopTree segmentloop_tree_05 = SegmentLoopTree();
	segmentloop_tree_05.insert(segmentloop_01);
	segmentloop_tree_05.insert(segmentloop_05);
	segmentloop_tree_05.insert(segmentloop_02);
	segmentloop_tree_05.insert(segmentloop_04);
	segmentloop_tree_05.insert(segmentloop_06);
	segmentloop_tree_05.insert(segmentloop_03);
	
	BOOST_CHECK(segmentloop_tree_05.search(p1) != nullptr);
	BOOST_CHECK(segmentloop_tree_05.search(p1)->gibbs_free_energy == -1.0);
	
	
	BOOST_CHECK(segmentloop_tree_05.size() == 6);
}

BOOST_AUTO_TEST_SUITE_END()
