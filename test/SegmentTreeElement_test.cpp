/**
 * @file test/SegmentTreeElement_test.cpp
 *
 * @date 2015-05-02
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



#define BOOST_TEST_MODULE SegmentTreeElement



#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "SubSequence.hpp"

#include "Segment.hpp"
#include "SegmentTreeElement.hpp"

#include <array>



#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests size of the SegmeentTreeElement
 *
 * @test
 *
 * @date 2015-05-01
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	std::string        segment_name  = "C/D-box K-turn";
	Sequence           sequence_5p   = Sequence("ACUUG");
	std::vector <Pair> bonds         = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p   = Sequence("AUG");
	float              energy        = -1.234;
	Segment            segment_01    = Segment(segment_name, sequence_5p, bonds, sequence_3p, energy);
	
	segment_name = "C/D-box K-turn";
	sequence_5p  = Sequence("ACUUG");
	bonds        = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	sequence_3p  = Sequence("AUC");
	Segment             segment_02   = Segment(segment_name, sequence_5p, bonds, sequence_3p, energy);
	
	SegmentTreeElement element_01 = SegmentTreeElement(segment_01);
	element_01.add_segment(segment_02);
	
	BOOST_CHECK_EQUAL(element_01.size(), 2);
}

/**
 * @brief Tests whether segments with a different 5' length can be searched through each other
 *
 * @test
 *
 * @date 2015-05-01
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	//ACUUG
	//:::
	//AUG
	std::string        segment_name_01 = "Segment 1";
	Sequence           sequence_5p_01  = Sequence("ACUUG");
	std::vector <Pair> bonds           = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_01  = Sequence("GUA");//Must be rotated to 5' -> 3'
	float              energy_01       = -1.00;
	Segment            segment_01      = Segment(segment_name_01, sequence_5p_01, bonds, sequence_3p_01, energy_01);
	Position           segment_01_5p_s = sequence_5p_01.data.begin();
	Position           segment_01_5p_e = sequence_5p_01.data.end() - 1;
	Position           segment_01_3p_s = sequence_3p_01.data.begin();
	Position           segment_01_3p_e = sequence_3p_01.data.end() - 1;
	SubSequence        segment_01_5p   = SubSequence(segment_01_5p_s, segment_01_5p_e);
	SubSequence        segment_01_3p   = SubSequence(segment_01_3p_s, segment_01_3p_e);
	
	//ACUU
	//:::
	//AUG
	std::string        segment_name_02 = "Segment 2";
	Sequence           sequence_5p_02  = Sequence("ACUU");
	bonds           = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_02  = Sequence("AUG");
	float              energy_02       = -2.00;
	Segment            segment_02      = Segment(segment_name_02, sequence_5p_02, bonds, sequence_3p_02, energy_02);
	SubSequence        segment_02_5p   = SubSequence(sequence_5p_02.data.begin(), (sequence_5p_02.data.end() - 1));
	SubSequence        segment_02_3p   = SubSequence(sequence_3p_02.data.begin(), (sequence_3p_02.data.end() - 1));
	
	// checks on initiazation
	BOOST_REQUIRE(sequence_5p_02.size() == 4);
	BOOST_REQUIRE(*(sequence_5p_02.data.begin() + 0) == Nucleotide::A);
	BOOST_REQUIRE(*(sequence_5p_02.data.begin() + 1) == Nucleotide::C);
	BOOST_REQUIRE(*(sequence_5p_02.data.begin() + 2) == Nucleotide::T);
	BOOST_REQUIRE(*(sequence_5p_02.data.begin() + 3) == Nucleotide::U);
	BOOST_REQUIRE((sequence_5p_02.data.begin() + 3) == (sequence_5p_02.data.end() - 1));
	
	BOOST_REQUIRE(sequence_3p_02.size() == 3);
	BOOST_REQUIRE(*(sequence_3p_02.data.begin() + 0) == Nucleotide::A);
	BOOST_REQUIRE(*(sequence_3p_02.data.begin() + 1) == Nucleotide::U);
	BOOST_REQUIRE(*(sequence_3p_02.data.begin() + 2) == Nucleotide::G);
	BOOST_REQUIRE((sequence_3p_02.data.begin() + 2) == (sequence_3p_02.data.end() - 1));
	
	
	SegmentTreeElement element_01 = SegmentTreeElement(segment_01);
	SegmentTreeElement element_02 = SegmentTreeElement(segment_02);
	element_01.add_segment(segment_02);
	element_02.add_segment(segment_01);
	
	Segment *segment_03 = element_01.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_04 = element_01.search_segment(segment_02_5p, segment_02_3p);
	Segment *segment_05 = element_02.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_06 = element_02.search_segment(segment_02_5p, segment_02_3p);
	
	// checks on the SegmentTreeElement
	BOOST_REQUIRE(segment_03 != nullptr);
	BOOST_REQUIRE(segment_04 != nullptr);
	BOOST_REQUIRE(segment_05 != nullptr);
	BOOST_REQUIRE(segment_06 != nullptr);
	
	BOOST_CHECK(segment_03->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_04->get_gibbs_free_energy() == energy_02);
	BOOST_CHECK(segment_05->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_06->get_gibbs_free_energy() == energy_02);
}

/**
 * @brief Tests whether segments with a different 5' sequence can be searched through each other
 *
 * @test
 *
 * @date 2015-05-01
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	//ACUUG
	//:::
	//AUG
	std::string        segment_name_01 = "Segment 1";
	Sequence           sequence_5p_01  = Sequence("ACUUG");
	std::vector <Pair> bonds           = {{Pair({0, 2}), Pair({2, 1}), Pair({4, 0})}};
	Sequence           sequence_3p_01  = Sequence("GUA");//rotate
	float              energy_01       = -1.00;
	Segment            segment_01      = Segment(segment_name_01, sequence_5p_01, bonds, sequence_3p_01, energy_01);
	Position           segment_01_5p_s = sequence_5p_01.data.begin();
	Position           segment_01_5p_e = sequence_5p_01.data.end() - 1;
	Position           segment_01_3p_s = sequence_3p_01.data.begin();
	Position           segment_01_3p_e = sequence_3p_01.data.end() - 1;
	SubSequence        segment_01_5p   = SubSequence(segment_01_5p_s, segment_01_5p_e);
	SubSequence        segment_01_3p   = SubSequence(segment_01_3p_s, segment_01_3p_e);
	
	//ACUUC
	//:::
	//AUG
	std::string        segment_name_02 = "Segment 2";
	Sequence           sequence_5p_02  = Sequence("ACUUC");
	bonds           = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0})}};
	Sequence           sequence_3p_02  = Sequence("GUA");
	float              energy_02       = -2.00;
	Segment            segment_02      = Segment(segment_name_02, sequence_5p_02, bonds, sequence_3p_02, energy_02);
	Position           segment_02_5p_s = sequence_5p_02.data.begin();
	Position           segment_02_5p_e = sequence_5p_02.data.end() - 1;
	Position           segment_02_3p_s = sequence_3p_02.data.begin();
	Position           segment_02_3p_e = sequence_3p_02.data.end() - 1;
	SubSequence        segment_02_5p   = SubSequence(segment_02_5p_s, segment_02_5p_e);
	SubSequence        segment_02_3p   = SubSequence(segment_02_3p_s, segment_02_3p_e);
	
	
	SegmentTreeElement element_01 = SegmentTreeElement(segment_01);
	SegmentTreeElement element_02 = SegmentTreeElement(segment_02);
	element_01.add_segment(segment_02);
	element_02.add_segment(segment_01);
	
	Segment *segment_03 = element_01.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_04 = element_01.search_segment(segment_02_5p, segment_02_3p);
	Segment *segment_05 = element_02.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_06 = element_02.search_segment(segment_02_5p, segment_02_3p);
	
	BOOST_REQUIRE(segment_03 != nullptr);
	BOOST_REQUIRE(segment_04 != nullptr);
	BOOST_REQUIRE(segment_05 != nullptr);
	BOOST_REQUIRE(segment_06 != nullptr);
	
	BOOST_CHECK(segment_03->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_04->get_gibbs_free_energy() == energy_02);
	BOOST_CHECK(segment_05->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_06->get_gibbs_free_energy() == energy_02);
}


/**
 * @brief Tests whether segments with a different 3' length can be searched through each other
 *
 * @test
 *
 * @date 2015-05-01
 */
BOOST_AUTO_TEST_CASE(Test4)
{
	//ACUUG
	//:::
	//AUG
	std::string        segment_name_01 = "Segment 1";
	Sequence           sequence_5p_01  = Sequence("ACUUG");
	std::vector <Pair> bonds           = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_01  = Sequence("GUA");// Rotate to 5'->3'
	float              energy_01       = -1.00;
	Segment            segment_01      = Segment(segment_name_01, sequence_5p_01, bonds, sequence_3p_01, energy_01);
	Position           segment_01_5p_s = sequence_5p_01.data.begin();
	Position           segment_01_5p_e = sequence_5p_01.data.end() - 1;
	Position           segment_01_3p_s = sequence_3p_01.data.begin();
	Position           segment_01_3p_e = sequence_3p_01.data.end() - 1;
	SubSequence        segment_01_5p   = SubSequence(segment_01_5p_s, segment_01_5p_e);
	SubSequence        segment_01_3p   = SubSequence(segment_01_3p_s, segment_01_3p_e);
	
	//ACUUG
	//:::
	//AUGA
	std::string        segment_name_02 = "Segment 2";
	Sequence           sequence_5p_02  = Sequence("ACUUG");
	bonds           = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_02  = Sequence("AGUA");// Rotate to 5'->3'
	float              energy_02       = -2.00;
	Segment            segment_02      = Segment(segment_name_02, sequence_5p_02, bonds, sequence_3p_02, energy_02);
	Position           segment_02_5p_s = sequence_5p_02.data.begin();
	Position           segment_02_5p_e = sequence_5p_02.data.end() - 1;
	Position           segment_02_3p_s = sequence_3p_02.data.begin();
	Position           segment_02_3p_e = sequence_3p_02.data.end() - 1;
	SubSequence        segment_02_5p   = SubSequence(segment_02_5p_s, segment_02_5p_e);
	SubSequence        segment_02_3p   = SubSequence(segment_02_3p_s, segment_02_3p_e);
	
	
	SegmentTreeElement element_01 = SegmentTreeElement(segment_01);
	SegmentTreeElement element_02 = SegmentTreeElement(segment_02);
	element_01.add_segment(segment_02);
	element_02.add_segment(segment_01);
	
	Segment *segment_03 = element_01.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_04 = element_01.search_segment(segment_02_5p, segment_02_3p);
	Segment *segment_05 = element_02.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_06 = element_02.search_segment(segment_02_5p, segment_02_3p);
	
	BOOST_REQUIRE(segment_03 != nullptr);
	BOOST_REQUIRE(segment_04 != nullptr);
	BOOST_REQUIRE(segment_05 != nullptr);
	BOOST_REQUIRE(segment_06 != nullptr);
	
	BOOST_CHECK(segment_03->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_04->get_gibbs_free_energy() == energy_02);
	BOOST_CHECK(segment_05->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_06->get_gibbs_free_energy() == energy_02);
}

/**
 * @brief Tests whether segments with a different 3' sequence can be searched through each other
 *
 * @test
 *
 * @date 2015-05-01
 */
BOOST_AUTO_TEST_CASE(Test5)
{
	//ACUUG
	//:::
	//AUG
	std::string        segment_name_01 = "Segment 1";
	Sequence           sequence_5p_01  = Sequence("ACUUG");
	std::vector <Pair> bonds_01        = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_01  = Sequence("GUA");// Reversing from 5' to 3' is essential
	float              energy_01       = -1.00;
	Segment            segment_01      = Segment(segment_name_01, sequence_5p_01, bonds_01, sequence_3p_01, energy_01);
	Position           segment_01_5p_s = sequence_5p_01.data.begin();
	Position           segment_01_5p_e = sequence_5p_01.data.end() - 1;
	Position           segment_01_3p_s = sequence_3p_01.data.begin();
	Position           segment_01_3p_e = sequence_3p_01.data.end() - 1;
	SubSequence        segment_01_5p   = SubSequence(segment_01_5p_s, segment_01_5p_e);
	SubSequence        segment_01_3p   = SubSequence(segment_01_3p_s, segment_01_3p_e);
	
	//ACUUG
	//:::
	//AUC
	std::string        segment_name_02 = "Segment 2";
	Sequence           sequence_5p_02  = Sequence("ACUUG");
	std::vector <Pair> bonds_02        = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_02  = Sequence("CUA");// Reversing from 5' to 3' is essential
	float              energy_02       = -2.00;
	Segment            segment_02      = Segment(segment_name_02, sequence_5p_02, bonds_02, sequence_3p_02, energy_02);
	Position           segment_02_5p_s = sequence_5p_02.data.begin();
	Position           segment_02_5p_e = sequence_5p_02.data.end() - 1;
	Position           segment_02_3p_s = sequence_3p_02.data.begin();
	Position           segment_02_3p_e = sequence_3p_02.data.end() - 1;
	SubSequence        segment_02_5p   = SubSequence(segment_02_5p_s, segment_02_5p_e);
	SubSequence        segment_02_3p   = SubSequence(segment_02_3p_s, segment_02_3p_e);
	
	
	SegmentTreeElement element_01 = SegmentTreeElement(segment_01);
	SegmentTreeElement element_02 = SegmentTreeElement(segment_02);
	element_01.add_segment(segment_02);
	element_02.add_segment(segment_01);
	
	Segment *segment_03 = element_01.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_04 = element_01.search_segment(segment_02_5p, segment_02_3p);
	Segment *segment_05 = element_02.search_segment(segment_01_5p, segment_01_3p);
	Segment *segment_06 = element_02.search_segment(segment_02_5p, segment_02_3p);
	
	BOOST_REQUIRE(segment_03 != nullptr);
	BOOST_REQUIRE(segment_04 != nullptr);
	BOOST_REQUIRE(segment_05 != nullptr);
	BOOST_REQUIRE(segment_06 != nullptr);
	
	BOOST_CHECK(segment_03->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_04->get_gibbs_free_energy() == energy_02);
	BOOST_CHECK(segment_05->get_gibbs_free_energy() == energy_01);
	BOOST_CHECK(segment_06->get_gibbs_free_energy() == energy_02);
}

/**
 * @brief Tests error throwing when inserting segments with identical sequences
 *
 * @test
 *
 * @date 2015-05-01
 */
BOOST_AUTO_TEST_CASE(Test6)
{
	std::string        segment_name_01 = "Segment 1";
	Sequence           sequence_5p_01  = Sequence("ACUUG");
	std::vector <Pair> bonds           = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_01  = Sequence("AUG");
	float              energy_01       = -1.00;
	Segment            segment_01      = Segment(segment_name_01, sequence_5p_01, bonds, sequence_3p_01, energy_01);
	
	std::string        segment_name_02 = "Segment 2";
	Sequence           sequence_5p_02  = Sequence("ACUUG");
	bonds           = {{ Pair({0, 2}), Pair({2, 1}), Pair({4, 0}) }};
	Sequence           sequence_3p_02  = Sequence("AUG");
	float              energy_02       = -2.00;
	Segment            segment_02      = Segment(segment_name_02, sequence_5p_02, bonds, sequence_3p_02, energy_02);
	
	BOOST_REQUIRE_EQUAL(sequence_5p_01.str(), sequence_5p_02.str());
	BOOST_REQUIRE_EQUAL(sequence_3p_01.str(), sequence_3p_02.str());
	
	BOOST_REQUIRE(sequence_5p_01.str() != sequence_3p_01.str());
	BOOST_REQUIRE(sequence_5p_01.str() != sequence_3p_02.str());
	BOOST_REQUIRE(sequence_3p_01.str() != sequence_5p_01.str());
	BOOST_REQUIRE(sequence_3p_01.str() != sequence_5p_02.str());
	
	SegmentTreeElement element_01 = SegmentTreeElement(segment_01);
	SegmentTreeElement element_02 = SegmentTreeElement(segment_02);
	
	BOOST_CHECK_THROW(element_01.add_segment(segment_02), std::invalid_argument);
	BOOST_CHECK_THROW(element_02.add_segment(segment_01), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
