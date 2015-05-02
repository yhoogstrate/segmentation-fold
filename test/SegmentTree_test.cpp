/**
 * @file test/SegmentTree_test.cpp
 *
 * @date 01-may-2015
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



#define BOOST_TEST_MODULE SegmentTree


#include "main.hpp"

#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Position.hpp"
#include "PairingPlus.hpp"
#include "Sequence.hpp"

#include "Segment.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"

#include <array>
#include <vector>

#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests quering using a PairingPlus element instead of deep copied sub sequences
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	SegmentTree segment_tree = SegmentTree();
	
	BOOST_CHECK(segment_tree.empty());
	BOOST_CHECK_EQUAL(segment_tree.size(), 0);
	
	
	std::string segment_name = "Articificial example";
	
	// ACUUG
	// \|/
	// AUG
	
	Sequence sequence_5p = Sequence("ACUUG");
	std::vector <Pair> bonds = {{Pair({2, 2}), Pair({4, 1}), Pair({5, 0})}};
	Sequence sequence_3p = Sequence("GUA");// rotation must be in the correct (5' -> 3') direction
	
	float energy = -1.234;
	
	Segment segment = Segment(segment_name, sequence_5p, bonds, sequence_3p, energy);
	
	
	segment_tree.insert(segment);
	
	BOOST_CHECK(!segment_tree.empty());
	BOOST_CHECK_EQUAL(segment_tree.size(), 1);
	
	
	
	// Setup query
	Sequence sequence = Sequence("ggg" "ACUUG" "ggg" "aaa" "ccc" "GUA" "ccc");
	
	Position p1a = sequence.data.begin() + 3;
	Position p1b = sequence.data.begin() + 7;
	
	BOOST_CHECK_EQUAL(*p1a , Nucleotide::A);
	BOOST_CHECK_EQUAL(*p1b , Nucleotide::G);
	
	
	Position p2a = sequence.data.begin() + 17;
	Position p2b = sequence.data.begin() + 19;
	
	BOOST_CHECK_EQUAL(*p2a , Nucleotide::G);
	BOOST_CHECK_EQUAL(*p2b , Nucleotide::A);
	
	
	PairingPlus p1 = PairingPlus(p1a, p1b);
	PairingPlus p2 = PairingPlus(p2a, p2b);
	
	
	Segment *segment_query = segment_tree.search(p1, p2);
	
	BOOST_CHECK(segment_query != nullptr);
}

/**
 * @brief Tests if searching a 0 size tree properly returns a nullptr
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	//boost::unit_test::unit_test_log.set_threshold_level( boost::unit_test::log_warnings );// enable warnings
	
	SegmentTree segment_tree = SegmentTree();
	
	BOOST_CHECK(segment_tree.empty());
	BOOST_CHECK_EQUAL(segment_tree.size(), 0);
	
	Sequence sequence = Sequence("cccACUUGcccaaagggGUAggg");
	
	
	
	Position p1a = sequence.data.begin() + 3;
	Position p1b = sequence.data.begin() + 7;
	
	BOOST_CHECK_EQUAL(*p1a , Nucleotide::A);
	BOOST_CHECK_EQUAL(*p1b , Nucleotide::G);
	
	
	Position p2a = sequence.data.begin() + 17;
	Position p2b = sequence.data.begin() + 19;
	
	BOOST_CHECK_EQUAL(*p2a , Nucleotide::G);
	BOOST_CHECK_EQUAL(*p2b , Nucleotide::A);
	
	
	PairingPlus subsequence_5p = PairingPlus(p1a, p1b);
	PairingPlus subsequence_3p = PairingPlus(p2a, p2b);
	
	
	Segment *segment_query01 = segment_tree.search(subsequence_5p, subsequence_3p);
	
	BOOST_CHECK(segment_query01 == nullptr);
}

/**
 * @brief Tests if for a tree with 1 element searching works
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	//boost::unit_test::unit_test_log.set_threshold_level(boost::unit_test::log_warnings);  // enable warnings
	
	std::string        segment_01_name  = "Segment 1";
	Sequence           segment_01_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_01_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0})}};
	Sequence           segment_01_seq3p = Sequence("GUA");
	float              segment_01_nrg   = -1.234;
	Segment            segment_01       = Segment(segment_01_name, segment_01_seq5p, segment_01_bonds, segment_01_seq3p, segment_01_nrg);
	
	SegmentTree segment_tree = SegmentTree();
	
	BOOST_CHECK(segment_tree.empty());
	BOOST_CHECK_EQUAL(segment_tree.size(), 0);
	
	segment_tree.insert(segment_01);
	
	BOOST_CHECK_EQUAL(segment_tree.size(), 1);
	
	
	
	Position p1a = segment_01_seq5p.data.begin();
	Position p1b = segment_01_seq5p.data.end() - 1;
	
	BOOST_CHECK_EQUAL(*p1a , Nucleotide::A);
	BOOST_CHECK_EQUAL(*p1b , Nucleotide::G);
	
	
	Position p2a = segment_01_seq3p.data.begin();
	Position p2b = segment_01_seq3p.data.end() - 1;
	
	
	BOOST_CHECK_EQUAL(*p2a , Nucleotide::G);
	BOOST_CHECK_EQUAL(*p2b , Nucleotide::A);
	
	
	PairingPlus subsequence_5p = PairingPlus(p1a, p1b);
	PairingPlus subsequence_3p = PairingPlus(p2a, p2b);
	
	Segment *segment_query01 = segment_tree.search(subsequence_5p, subsequence_3p);
	
	BOOST_CHECK(segment_query01 != nullptr);
	BOOST_CHECK(segment_01.name == segment_query01->name);
	
	// See if the original object is returned instead of a (new) deep copy
	BOOST_CHECK(&segment_01 == segment_query01);
}

/**
 * @brief Tests if 2 sized tree searching works
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test4)
{
	boost::unit_test::unit_test_log.set_threshold_level(boost::unit_test::log_warnings);  // enable warnings
	
	std::string        segment_01_name  = "Segment 1";
	Sequence           segment_01_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_01_bonds = {{Pair({0, 2}), Pair({2, 1}), Pair({4, 0})}};
	Sequence           segment_01_seq3p = Sequence("AUG");
	float              segment_01_nrg   = -1.234;
	Segment            segment_01       = Segment(segment_01_name, segment_01_seq5p, segment_01_bonds, segment_01_seq3p, segment_01_nrg);
	
	std::string        segment_02_name  = "Segment 2";
	Sequence           segment_02_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_02_bonds = {{Pair({0, 2}), Pair({2, 1}), Pair({4, 0})}};
	Sequence           segment_02_seq3p = Sequence("AUC");
	float              segment_02_nrg   = -1.234;
	Segment            segment_02       = Segment(segment_02_name, segment_02_seq5p, segment_02_bonds, segment_02_seq3p, segment_02_nrg);
	
	SegmentTree segment_tree = SegmentTree();
	
	BOOST_CHECK(segment_tree.empty() == true);
	BOOST_CHECK(segment_tree.size() == 0);
	
	// Add 1
	segment_tree.insert(segment_01);
	BOOST_CHECK(segment_tree.empty() == false);
	BOOST_CHECK(segment_tree.size() == 1);
	
	// Add 2
	segment_tree.insert(segment_02);
	BOOST_CHECK(segment_tree.empty() == false);
	BOOST_CHECK(segment_tree.size() == 2);
	
	// Obtain 1
	Position p1a = segment_01_seq5p.data.begin();
	Position p1b = segment_01_seq5p.data.end() - 1;
	Position p2a = segment_01_seq3p.data.begin();
	Position p2b = segment_01_seq3p.data.end() - 1;
	PairingPlus subsequence1_5p = PairingPlus(p1a, p1b);
	PairingPlus subsequence1_3p = PairingPlus(p2a, p2b);
	
	Segment *segment_query01 = segment_tree.search(subsequence1_5p, subsequence1_3p);
	BOOST_CHECK(segment_query01 != nullptr);
	BOOST_CHECK(segment_01.name == segment_query01->name);
	BOOST_CHECK(&segment_01 == segment_query01);// Deep copy detection
	
	// Obtain 2
	Position q1a = segment_02_seq5p.data.begin();
	Position q1b = segment_02_seq5p.data.end() - 1;
	Position q2a = segment_02_seq3p.data.begin();
	Position q2b = segment_02_seq3p.data.end() - 1;
	PairingPlus subsequence2_5p = PairingPlus(q1a, q1b);
	PairingPlus subsequence2_3p = PairingPlus(q2a, q2b);
	
	Segment *segment_query02 = segment_tree.search(subsequence2_5p, subsequence2_3p);
	BOOST_CHECK(segment_query02 != nullptr);
	BOOST_CHECK(segment_02.name == segment_query02->name);
	BOOST_CHECK(&segment_02 == segment_query02);// Deep copy detection
	
	// Further checks
	BOOST_CHECK(segment_01.name != segment_query02->name);
	BOOST_CHECK(segment_02.name != segment_query01->name);
}

/**
 * @brief Tests if for 5 sized trees (with elements of different sizes and equences) searching works, in all possible orders of being added
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test5)
{
	//ACUUG
	//:::
	//AUG
	std::string        segment_01_name  = "Segment 1";
	Sequence           segment_01_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_01_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_01_seq3p = Sequence("GUA");// must be given in 5' to 3' orientation
	float              segment_01_nrg   = -1.0;
	Segment            segment_01       = Segment(segment_01_name, segment_01_seq5p, segment_01_bonds, segment_01_seq3p, segment_01_nrg);
	Position           segment_01_5p_s  = segment_01_seq5p.data.begin();
	Position           segment_01_5p_e  = segment_01_seq5p.data.end() - 1;
	Position           segment_01_3p_s  = segment_01_seq3p.data.begin();
	Position           segment_01_3p_e  = segment_01_seq3p.data.end() - 1;
	PairingPlus        segment_01_5p    = PairingPlus(segment_01_5p_s, segment_01_5p_e);
	PairingPlus        segment_01_3p    = PairingPlus(segment_01_3p_s, segment_01_3p_e);
	
	//ACUUG
	//:::
	//AUC
	std::string        segment_02_name  = "Segment 2";
	Sequence           segment_02_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_02_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_02_seq3p = Sequence("CUA");// must be given in 5' to 3' orientation
	float              segment_02_nrg   = -2.0;
	Segment            segment_02       = Segment(segment_02_name, segment_02_seq5p, segment_02_bonds, segment_02_seq3p, segment_02_nrg);
	Position           segment_02_5p_s  = segment_02_seq5p.data.begin();
	Position           segment_02_5p_e  = segment_02_seq5p.data.end() - 1;
	Position           segment_02_3p_s  = segment_02_seq3p.data.begin();
	Position           segment_02_3p_e  = segment_02_seq3p.data.end() - 1;
	PairingPlus        segment_02_5p    = PairingPlus(segment_02_5p_s, segment_02_5p_e);
	PairingPlus        segment_02_3p    = PairingPlus(segment_02_3p_s, segment_02_3p_e);
	
	//ACUU
	//:::
	//AUG
	std::string        segment_03_name  = "Segment 3";
	Sequence           segment_03_seq5p = Sequence("ACUU");
	std::vector <Pair> segment_03_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_03_seq3p = Sequence("GUA");// must be given in 5' to 3' orientation
	float              segment_03_nrg   = -3.0;
	Segment            segment_03       = Segment(segment_03_name, segment_03_seq5p, segment_03_bonds, segment_03_seq3p, segment_03_nrg);
	Position           segment_03_5p_s  = segment_03_seq5p.data.begin();
	Position           segment_03_5p_e  = segment_03_seq5p.data.end() - 1;
	Position           segment_03_3p_s  = segment_03_seq3p.data.begin();
	Position           segment_03_3p_e  = segment_03_seq3p.data.end() - 1;
	PairingPlus        segment_03_5p    = PairingPlus(segment_03_5p_s, segment_03_5p_e);
	PairingPlus        segment_03_3p    = PairingPlus(segment_03_3p_s, segment_03_3p_e);
	
	//ACUUC
	//:::
	//AUG
	std::string        segment_04_name  = "Segment 4";
	Sequence           segment_04_seq5p = Sequence("ACUUC");
	std::vector <Pair> segment_04_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_04_seq3p = Sequence("GUA");// must be given in 5' to 3' orientation
	float              segment_04_nrg   = -4.0;
	Segment            segment_04       = Segment(segment_04_name, segment_04_seq5p, segment_04_bonds, segment_04_seq3p, segment_04_nrg);
	Position           segment_04_5p_s  = segment_04_seq5p.data.begin();
	Position           segment_04_5p_e  = segment_04_seq5p.data.end() - 1;
	Position           segment_04_3p_s  = segment_04_seq3p.data.begin();
	Position           segment_04_3p_e  = segment_04_seq3p.data.end() - 1;
	PairingPlus        segment_04_5p    = PairingPlus(segment_04_5p_s, segment_04_5p_e);
	PairingPlus        segment_04_3p    = PairingPlus(segment_04_3p_s, segment_04_3p_e);
	
	//ACUUG
	//:::
	//AUGA
	std::string        segment_05_name  = "Segment 5";
	Sequence           segment_05_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_05_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_05_seq3p = Sequence("AGUA");// must be given in 5' to 3' orientation
	float              segment_05_nrg   = -5.0;
	Segment            segment_05       = Segment(segment_05_name, segment_05_seq5p, segment_05_bonds, segment_05_seq3p, segment_05_nrg);
	Position           segment_05_5p_s  = segment_05_seq5p.data.begin();
	Position           segment_05_5p_e  = segment_05_seq5p.data.end() - 1;
	Position           segment_05_3p_s  = segment_05_seq3p.data.begin();
	Position           segment_05_3p_e  = segment_05_seq3p.data.end() - 1;
	PairingPlus        segment_05_5p    = PairingPlus(segment_05_5p_s, segment_05_5p_e);
	PairingPlus        segment_05_3p    = PairingPlus(segment_05_3p_s, segment_05_3p_e);
	
	char permutations[120][5] = {{1, 2, 3, 4, 5}, {1, 2, 3, 5, 4}, {1, 2, 4, 3, 5}, {1, 2, 4, 5, 3}, {1, 2, 5, 3, 4}, {1, 2, 5, 4, 3}, {1, 3, 2, 4, 5}, {1, 3, 2, 5, 4}, {1, 3, 4, 2, 5}, {1, 3, 4, 5, 2}, {1, 3, 5, 2, 4}, {1, 3, 5, 4, 2}, {1, 4, 2, 3, 5}, {1, 4, 2, 5, 3}, {1, 4, 3, 2, 5}, {1, 4, 3, 5, 2}, {1, 4, 5, 2, 3}, {1, 4, 5, 3, 2}, {1, 5, 2, 3, 4}, {1, 5, 2, 4, 3}, {1, 5, 3, 2, 4}, {1, 5, 3, 4, 2}, {1, 5, 4, 2, 3}, {1, 5, 4, 3, 2}, {2, 1, 3, 4, 5}, {2, 1, 3, 5, 4}, {2, 1, 4, 3, 5}, {2, 1, 4, 5, 3}, {2, 1, 5, 3, 4}, {2, 1, 5, 4, 3}, {2, 3, 1, 4, 5}, {2, 3, 1, 5, 4}, {2, 3, 4, 1, 5}, {2, 3, 4, 5, 1}, {2, 3, 5, 1, 4}, {2, 3, 5, 4, 1}, {2, 4, 1, 3, 5}, {2, 4, 1, 5, 3}, {2, 4, 3, 1, 5}, {2, 4, 3, 5, 1}, {2, 4, 5, 1, 3}, {2, 4, 5, 3, 1}, {2, 5, 1, 3, 4}, {2, 5, 1, 4, 3}, {2, 5, 3, 1, 4}, {2, 5, 3, 4, 1}, {2, 5, 4, 1, 3}, {2, 5, 4, 3, 1}, {3, 1, 2, 4, 5}, {3, 1, 2, 5, 4}, {3, 1, 4, 2, 5}, {3, 1, 4, 5, 2}, {3, 1, 5, 2, 4}, {3, 1, 5, 4, 2}, {3, 2, 1, 4, 5}, {3, 2, 1, 5, 4}, {3, 2, 4, 1, 5}, {3, 2, 4, 5, 1}, {3, 2, 5, 1, 4}, {3, 2, 5, 4, 1}, {3, 4, 1, 2, 5}, {3, 4, 1, 5, 2}, {3, 4, 2, 1, 5}, {3, 4, 2, 5, 1}, {3, 4, 5, 1, 2}, {3, 4, 5, 2, 1}, {3, 5, 1, 2, 4}, {3, 5, 1, 4, 2}, {3, 5, 2, 1, 4}, {3, 5, 2, 4, 1}, {3, 5, 4, 1, 2}, {3, 5, 4, 2, 1}, {4, 1, 2, 3, 5}, {4, 1, 2, 5, 3}, {4, 1, 3, 2, 5}, {4, 1, 3, 5, 2}, {4, 1, 5, 2, 3}, {4, 1, 5, 3, 2}, {4, 2, 1, 3, 5}, {4, 2, 1, 5, 3}, {4, 2, 3, 1, 5}, {4, 2, 3, 5, 1}, {4, 2, 5, 1, 3}, {4, 2, 5, 3, 1}, {4, 3, 1, 2, 5}, {4, 3, 1, 5, 2}, {4, 3, 2, 1, 5}, {4, 3, 2, 5, 1}, {4, 3, 5, 1, 2}, {4, 3, 5, 2, 1}, {4, 5, 1, 2, 3}, {4, 5, 1, 3, 2}, {4, 5, 2, 1, 3}, {4, 5, 2, 3, 1}, {4, 5, 3, 1, 2}, {4, 5, 3, 2, 1}, {5, 1, 2, 3, 4}, {5, 1, 2, 4, 3}, {5, 1, 3, 2, 4}, {5, 1, 3, 4, 2}, {5, 1, 4, 2, 3}, {5, 1, 4, 3, 2}, {5, 2, 1, 3, 4}, {5, 2, 1, 4, 3}, {5, 2, 3, 1, 4}, {5, 2, 3, 4, 1}, {5, 2, 4, 1, 3}, {5, 2, 4, 3, 1}, {5, 3, 1, 2, 4}, {5, 3, 1, 4, 2}, {5, 3, 2, 1, 4}, {5, 3, 2, 4, 1}, {5, 3, 4, 1, 2}, {5, 3, 4, 2, 1}, {5, 4, 1, 2, 3}, {5, 4, 1, 3, 2}, {5, 4, 2, 1, 3}, {5, 4, 2, 3, 1}, {5, 4, 3, 1, 2}, {5, 4, 3, 2, 1}};
	
	unsigned int i, j;
	
	SegmentTree segment_tree = SegmentTree();
	
	for(i = 0; i < 120; i++)											// Test with 120 different orders of adding segments
	{
		SegmentTree segment_tree = SegmentTree();
		
		BOOST_CHECK(segment_tree.empty() == true);
		BOOST_CHECK(segment_tree.size() == 0);
		
		for(j = 0; j < 5; j++)
		{
			switch(permutations[i][j])
			{
				case 1:
					segment_tree.insert(segment_01);
					break;
				case 2:
					segment_tree.insert(segment_02);
					break;
				case 3:
					segment_tree.insert(segment_03);
					break;
				case 4:
					segment_tree.insert(segment_04);
					break;
				case 5:
					segment_tree.insert(segment_05);
					break;
			}
			
			BOOST_CHECK(segment_tree.size() == (j + 1));
		}
		
		// Obtain 01
		Segment *segment_query01 = segment_tree.search(segment_01_5p , segment_01_3p);
		BOOST_REQUIRE_MESSAGE(segment_query01 != nullptr, "Failure obtaining segment 01 at iteration " << i);
		BOOST_CHECK(segment_01.name == segment_query01->name);
		BOOST_CHECK(&segment_01 == segment_query01);// Deep copy detection
		
		// Obtain 02
		Segment *segment_query02 = segment_tree.search(segment_02_5p , segment_02_3p);
		BOOST_REQUIRE_MESSAGE(segment_query02 != nullptr, "Failure obtaining segment 02 at iteration " << i);
		BOOST_CHECK(segment_02.name == segment_query02->name);
		BOOST_CHECK(&segment_02 == segment_query02);// Deep copy detection
		
		// Obtain 03
		Segment *segment_query03 = segment_tree.search(segment_03_5p , segment_03_3p);
		BOOST_REQUIRE_MESSAGE(segment_query03 != nullptr, "Failure obtaining segment 03 at iteration " << i);
		BOOST_CHECK(segment_03.name == segment_query03->name);
		BOOST_CHECK(&segment_03 == segment_query03);// Deep copy detection
		
		// Obtain 04
		Segment *segment_query04 = segment_tree.search(segment_04_5p , segment_04_3p);
		BOOST_REQUIRE_MESSAGE(segment_query04 != nullptr, "Failure obtaining segment 04 at iteration " << i);
		BOOST_CHECK(segment_04.name == segment_query04->name);
		BOOST_CHECK(&segment_04 == segment_query04);// Deep copy detection
		
		// Obtain 05
		Segment *segment_query05 = segment_tree.search(segment_05_5p , segment_05_3p);
		BOOST_REQUIRE_MESSAGE(segment_query05 != nullptr, "Failure obtaining segment 05 at iteration " << i);
		BOOST_CHECK(segment_05.name == segment_query05->name);
		BOOST_CHECK(&segment_05 == segment_query05);// Deep copy detection
	}
	
	// Request a non-existing sequence and checks whether that returns a nullptr
	Sequence           segment_xx_seq   = Sequence("AAAAAgggCCC");
	Position           segment_xx_5p_s  = segment_xx_seq.data.begin();
	Position           segment_xx_5p_e  = segment_xx_seq.data.begin() + 4;
	Position           segment_xx_3p_s  = segment_xx_seq.data.begin() + 4 + 3 + 1;
	Position           segment_xx_3p_e  = segment_xx_seq.data.begin() + 4 + 3 + 3;
	PairingPlus        segment_xx_5p    = PairingPlus(segment_xx_5p_s, segment_xx_5p_e); //AAAAA
	PairingPlus        segment_xx_3p    = PairingPlus(segment_xx_3p_s, segment_xx_3p_e); //CCC
	
	BOOST_CHECK(segment_tree.search(segment_xx_5p , segment_xx_3p) == nullptr);
}

/**
 * @brief Tests whether an exception is raised if a segment with a duplicate sequence is added
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test6)
{
	std::string        segment_01_name  = "Segment 1";
	Sequence           segment_01_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_01_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_01_seq3p = Sequence("AUG");
	float              segment_01_nrg   = -1.0;
	Segment            segment_01       = Segment(segment_01_name, segment_01_seq5p, segment_01_bonds, segment_01_seq3p, segment_01_nrg);
	
	std::string        segment_02_name  = "Segment 2";
	Sequence           segment_02_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_02_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_02_seq3p = Sequence("AUG");
	float              segment_02_nrg   = -2.0;
	Segment            segment_02       = Segment(segment_02_name, segment_02_seq5p, segment_02_bonds, segment_02_seq3p, segment_02_nrg);
	
	SegmentTree segment_tree = SegmentTree();
	
	segment_tree.insert(segment_01);
	BOOST_CHECK_THROW(segment_tree.insert(segment_02), std::invalid_argument);
}

/**
 * @brief Tests whether a reverse complement Segment does not throw an exception
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test7)
{
	std::string        segment_fwd_name  = "Segment 1 fwd";
	Sequence           segment_fwd_seq5p = Sequence("ACUUG");
	std::vector <Pair> segment_fwd_bonds = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           segment_fwd_seq3p = Sequence("AUG");
	float              segment_fwd_nrg   = -2.0;
	Segment            segment_fwd       = Segment(segment_fwd_name, segment_fwd_seq5p, segment_fwd_bonds, segment_fwd_seq3p, segment_fwd_nrg);
	Position           segment_fwd_5p_s  = segment_fwd_seq5p.data.begin();
	Position           segment_fwd_5p_e  = segment_fwd_seq5p.data.end() - 1;
	Position           segment_fwd_3p_s  = segment_fwd_seq3p.data.begin();
	Position           segment_fwd_3p_e  = segment_fwd_seq3p.data.end() - 1;
	PairingPlus        segment_fwd_5p    = PairingPlus(segment_fwd_5p_s, segment_fwd_5p_e);
	PairingPlus        segment_fwd_3p    = PairingPlus(segment_fwd_3p_s, segment_fwd_3p_e);
	
	std::string        segment_rev_name  = "Segment 1 rev";
	Sequence           segment_rev_seq5p = Sequence("GUA");
	std::vector <Pair> segment_rev_bonds = {{ Pair({0, 0}), Pair({1, 1}), Pair({2, 2}) }};
	Sequence           segment_rev_seq3p = Sequence("GUUCA");
	float              segment_rev_nrg   = -2.0;
	Segment            segment_rev       = Segment(segment_rev_name, segment_rev_seq5p, segment_rev_bonds, segment_rev_seq3p, segment_rev_nrg);
	Position           segment_rev_5p_s  = segment_rev_seq5p.data.begin();
	Position           segment_rev_5p_e  = segment_rev_seq5p.data.end() - 1;
	Position           segment_rev_3p_s  = segment_rev_seq3p.data.begin();
	Position           segment_rev_3p_e  = segment_rev_seq3p.data.end() - 1;
	PairingPlus        segment_rev_5p    = PairingPlus(segment_rev_5p_s, segment_rev_5p_e);
	PairingPlus        segment_rev_3p    = PairingPlus(segment_rev_3p_s, segment_rev_3p_e);
	
	SegmentTree segment_tree = SegmentTree();
	segment_tree.insert(segment_fwd);
	segment_tree.insert(segment_rev);
	
	BOOST_CHECK(segment_tree.search(segment_fwd_5p , segment_fwd_3p) != nullptr);
	BOOST_CHECK(segment_tree.search(segment_rev_5p , segment_rev_3p) != nullptr);
	BOOST_CHECK(segment_tree.search(segment_fwd_5p , segment_fwd_3p) != segment_tree.search(segment_rev_5p , segment_rev_3p));
	
	BOOST_CHECK_THROW(segment_tree.insert(segment_fwd), std::invalid_argument);
	BOOST_CHECK_THROW(segment_tree.insert(segment_rev), std::invalid_argument);
}
BOOST_AUTO_TEST_SUITE_END()
