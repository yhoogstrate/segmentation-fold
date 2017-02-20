/**
 * @file test/Sequence_test.cpp
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



#define BOOST_TEST_MODULE Sequence



#include "Nucleotide.hpp"
#include "Sequence.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether single nucleotides are inserted correctly, and whether the size of the sequence fits.
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test01)
{
	Sequence s;
	
	s = Sequence("ACTG");
	BOOST_CHECK_EQUAL(s.size(), 4);
	
	Nucleotide n1 = Nucleotide::A;
	const char n2 = 'C';
	char n3 = 'u';
	
	s.push_back(n1);
	s.push_back(n2);
	s.push_back(n3);
	s.push_back('G');
	
	BOOST_CHECK_EQUAL(s.size(), 8);
	
	BOOST_CHECK_EQUAL(s[0], Nucleotide::A);
	BOOST_CHECK_EQUAL(s[1], Nucleotide::C);
	BOOST_CHECK_EQUAL(s[2], Nucleotide::U);
	BOOST_CHECK_EQUAL(s[3], Nucleotide::G);
	
	BOOST_CHECK_EQUAL(s[4], Nucleotide::A);
	BOOST_CHECK_EQUAL(s[5], Nucleotide::C);
	BOOST_CHECK_EQUAL(s[6], Nucleotide::U);
	BOOST_CHECK_EQUAL(s[7], Nucleotide::G);
}



/**
 * @brief Checks std::string to Sequence and Sequence to std::string conversion
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test02)
{
	std::string s2_str = std::string("ACUG");
	Sequence s1 = Sequence("ACGU");
	Sequence s2 = Sequence(s2_str);
	
	Sequence s3 = Sequence(s1);
	Sequence s4 = Sequence(s2);
	
	BOOST_CHECK_EQUAL(s1.size(), 4);
	BOOST_CHECK_EQUAL(s2.size(), 4);
	BOOST_CHECK_EQUAL(s3.size(), 4);
	BOOST_CHECK_EQUAL(s4.size(), 4);
	
	BOOST_CHECK(s3.str().compare(s4.str()) != 0);
	BOOST_CHECK(s2_str.compare(s4.str()) == 0);
}



/**
 * @brief Tests behaviour of unexpected bases char-wise: should throw an exception
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test03)
{
	Sequence sequence;
	signed char c;
	
	for(c = -128; c < 127; c++)
	{
		sequence = Sequence();
		
		switch(c)
		{
			case 'a':
			case 'A':
			case 'c':
			case 'C':
			case 't':
			case 'T':
			case 'u':
			case 'U':
			case 'g':
			case 'G':
				BOOST_CHECK_NO_THROW(sequence.push_back(c));
				break;
			default:
				BOOST_CHECK_THROW(sequence.push_back(c), std::invalid_argument);
				break;
		}
	}
}



/**
 * @brief Tests behaviour of unexpected bases using entire strings: should throw an exception
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test04)
{
	Sequence sequence;
	std::string sequence_str;
	signed char c;
	
	for(c = -128; c < 127; c++)
	{
		sequence_str = "ActuG";
		sequence_str += c;
		
		switch(c)
		{
			case 'a':
			case 'A':
			case 'c':
			case 'C':
			case 't':
			case 'T':
			case 'u':
			case 'U':
			case 'g':
			case 'G':
				BOOST_CHECK_NO_THROW(sequence = Sequence(sequence_str));
				break;
			default:
				BOOST_CHECK_THROW(sequence = Sequence(sequence_str), std::invalid_argument);
				break;
		}
	}
}



/**
 *
 * @test
 *
 * @section DESCRIPTION
 * Given that nucleotides are enumerated as follows:
 * A = 0
 * C = 1
 * G = 2
 * U = 3
 * A sequence of AA must be the 'smallest' of the two-letter sequences.
 */
BOOST_AUTO_TEST_CASE(Test05)
{
	Sequence s1  = Sequence("AA");
	Sequence s2  = Sequence("CC");
	Sequence s3  = Sequence("CG");
	Sequence s4a = Sequence("GU");
	Sequence s4b = Sequence("GT");
	Sequence s5  = Sequence("TU");
	
	
	BOOST_CHECK(s1  < s2);
	BOOST_CHECK(s2  < s3);
	BOOST_CHECK(s3  < s4a);
	BOOST_CHECK(s3  < s4b);
	BOOST_CHECK(s4a < s5);
	BOOST_CHECK(s4b < s5);
	
	BOOST_CHECK(s1  < s3);
	BOOST_CHECK(s2  < s4a);
	BOOST_CHECK(s2  < s4b);
	BOOST_CHECK(s3  < s5);
	
	BOOST_CHECK(s1  < s4a);
	BOOST_CHECK(s1  < s4b);
	BOOST_CHECK(s2  < s5);
	
	BOOST_CHECK(s1  < s5);
	
	
	BOOST_CHECK((s2  < s1) == false);
	BOOST_CHECK((s3  < s2) == false);
	BOOST_CHECK((s4a < s3) == false);
	BOOST_CHECK((s4b < s3) == false);
	BOOST_CHECK((s5  < s4a) == false);
	BOOST_CHECK((s5  < s4b) == false);
	
	BOOST_CHECK((s3  < s1) == false);
	BOOST_CHECK((s4a < s2) == false);
	BOOST_CHECK((s4b < s2) == false);
	BOOST_CHECK((s5  < s3) == false);
	
	BOOST_CHECK((s4a < s1) == false);
	BOOST_CHECK((s4b < s1) == false);
	BOOST_CHECK((s5  < s2) == false);
	
	BOOST_CHECK((s5  < s1) == false);
}



/**
 * @brief Tests the '==' operator
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test06)
{
	Sequence s1 = Sequence("AAGAA");
	Sequence s2 = Sequence("aAGAA");
	Sequence s3 = Sequence("AAGA");
	Sequence s4 = Sequence("AAGAAA");
	Sequence s5 = Sequence("TTGTT");
	Sequence s6 = Sequence("TUGuu");
	
	
	BOOST_CHECK(s1  == s2);
	BOOST_CHECK((s2 == s3) == false);
	BOOST_CHECK((s3 == s4) == false);
	BOOST_CHECK((s4 == s5) == false);
	BOOST_CHECK(s5  == s6);
	
	BOOST_CHECK((s1 == s3) == false);
	BOOST_CHECK((s2 == s4) == false);
	BOOST_CHECK((s3 == s5) == false);
	BOOST_CHECK((s4 == s6) == false);
	
	BOOST_CHECK((s1 == s4) == false);
	BOOST_CHECK((s2 == s5) == false);
	BOOST_CHECK((s3 == s6) == false);
	
	BOOST_CHECK((s1 == s5) == false);
	BOOST_CHECK((s2 == s6) == false);
	
	BOOST_CHECK((s1 == s6) == false);
}



/**
 * @brief Tests the subseq function
 */
BOOST_AUTO_TEST_CASE(Test07)
{
	Sequence sequence = Sequence("ACTG");
	
	BOOST_REQUIRE_EQUAL(sequence.size(), 4);
	
	Sequence subsequence_01 = sequence.subseq(0, 0);// start at 0, stop at 0 -> A
	Sequence subsequence_02 = sequence.subseq(0, 1);// start at 0, stop at 1 -> AC
	Sequence subsequence_03 = sequence.subseq(0, 2);// start at 0, stop at 2 -> ACT
	Sequence subsequence_04 = sequence.subseq(0, 3);// start at 0, stop at 3 -> ACTG
	Sequence subsequence_05 = sequence.subseq(1, 3);// start at 1, stop at 3 ->  CTG
	Sequence subsequence_06 = sequence.subseq(2, 3);// start at 2, stop at 3 ->   TG
	Sequence subsequence_07 = sequence.subseq(3, 3);// start at 3, stop at 3 ->    G
	
	BOOST_CHECK_EQUAL(subsequence_01.size(), 1);
	BOOST_CHECK_EQUAL(subsequence_02.size(), 2);
	BOOST_CHECK_EQUAL(subsequence_03.size(), 3);
	BOOST_CHECK_EQUAL(subsequence_04.size(), 4);
	BOOST_CHECK_EQUAL(subsequence_05.size(), 3);
	BOOST_CHECK_EQUAL(subsequence_06.size(), 2);
	BOOST_CHECK_EQUAL(subsequence_07.size(), 1);
	
	BOOST_CHECK(subsequence_01 == Sequence("A"));
	BOOST_CHECK(subsequence_02 == Sequence("AC"));
	BOOST_CHECK(subsequence_03 == Sequence("ACT"));
	BOOST_CHECK(subsequence_04 == Sequence("ACTG"));
	BOOST_CHECK(subsequence_05 == Sequence("CTG"));
	BOOST_CHECK(subsequence_06 == Sequence("TG"));
	BOOST_CHECK(subsequence_07 == Sequence("G"));
	
	BOOST_CHECK((subsequence_01 == Sequence("AA")) == false);
	BOOST_CHECK((subsequence_02 == Sequence("ACA")) == false);
	BOOST_CHECK((subsequence_03 == Sequence("ACTA")) == false);
	BOOST_CHECK((subsequence_04 == Sequence("ACTGA")) == false);
	BOOST_CHECK((subsequence_05 == Sequence("CTGA")) == false);
	BOOST_CHECK((subsequence_06 == Sequence("TGA")) == false);
	BOOST_CHECK((subsequence_07 == Sequence("GA")) == false);
}



/**
 * @brief Valgrind check
 */
BOOST_AUTO_TEST_CASE(Test08)
{
	Sequence sequence = Sequence();
	BOOST_REQUIRE_EQUAL(sequence.size(), 0);
}



/**
 * @brief Tests Sequence::compare(&Sequence)
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test09)
{
	/*
	A = 0,
	C = 1,
	G = 2,
	U, T = 3
	*/
	
	Sequence sequence = Sequence("CCC");
	
	Sequence smaller_sequence_01 = Sequence("ACC");
	Sequence smaller_sequence_02 = Sequence("CCA");
	Sequence smaller_sequence_03 = Sequence("AA");
	Sequence smaller_sequence_04 = Sequence("CC");
	Sequence smaller_sequence_05 = Sequence("TT");
	
	Sequence larger_sequence_01 = Sequence("ACCC");
	Sequence larger_sequence_02 = Sequence("CCCA");
	Sequence larger_sequence_03 = Sequence("CCG");
	Sequence larger_sequence_04 = Sequence("GCC");
	Sequence larger_sequence_05 = Sequence("TTT");
	
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_sequence_01), IS_LARGER);// left sequence is larger than the 'smaller_' sequence
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_sequence_02), IS_LARGER);
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_sequence_03), IS_LARGER);
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_sequence_04), IS_LARGER);
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_sequence_05), IS_LARGER);
	
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_sequence_01), IS_SMALLER);// left sequence is larger than the 'larger_' sequence
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_sequence_02), IS_SMALLER);
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_sequence_03), IS_SMALLER);
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_sequence_04), IS_SMALLER);
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_sequence_05), IS_SMALLER);
}



/**
 * @brief Tests Sequence::ssubseq and Sequence::compare(&SubSequence)
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test10)
{
	Sequence sequence = Sequence("CCC");
	
	Sequence smaller_sequence_01 = Sequence("ACC");
	Sequence smaller_sequence_02 = Sequence("CCA");
	Sequence smaller_sequence_03 = Sequence("AA");
	Sequence smaller_sequence_04 = Sequence("CC");
	Sequence smaller_sequence_05 = Sequence("TT");
	
	SubSequence smaller_subsequence_01 = smaller_sequence_01.ssubseq(0, smaller_sequence_01.size() - 1);
	SubSequence smaller_subsequence_02 = smaller_sequence_02.ssubseq(0, smaller_sequence_02.size() - 1);
	SubSequence smaller_subsequence_03 = smaller_sequence_03.ssubseq(0, smaller_sequence_03.size() - 1);
	SubSequence smaller_subsequence_04 = smaller_sequence_04.ssubseq(0, smaller_sequence_04.size() - 1);
	SubSequence smaller_subsequence_05 = smaller_sequence_05.ssubseq(0, smaller_sequence_05.size() - 1);
	
	BOOST_CHECK(smaller_subsequence_01[0] == Nucleotide::A);
	BOOST_CHECK(smaller_subsequence_01[1] == Nucleotide::C);
	BOOST_CHECK(smaller_subsequence_01[2] == Nucleotide::C);
	BOOST_CHECK(smaller_subsequence_01.size == 3);
	
	
	Sequence larger_sequence_01 = Sequence("ACCC");
	Sequence larger_sequence_02 = Sequence("CCCA");
	Sequence larger_sequence_03 = Sequence("CCG");
	Sequence larger_sequence_04 = Sequence("GCC");
	Sequence larger_sequence_05 = Sequence("TTT");
	
	SubSequence larger_subsequence_01 = larger_sequence_01.ssubseq(0, larger_sequence_01.size() - 1);
	SubSequence larger_subsequence_02 = larger_sequence_02.ssubseq(0, larger_sequence_02.size() - 1);
	SubSequence larger_subsequence_03 = larger_sequence_03.ssubseq(0, larger_sequence_03.size() - 1);
	SubSequence larger_subsequence_04 = larger_sequence_04.ssubseq(0, larger_sequence_04.size() - 1);
	SubSequence larger_subsequence_05 = larger_sequence_05.ssubseq(0, larger_sequence_05.size() - 1);
	
	BOOST_CHECK(larger_subsequence_01[0] == Nucleotide::A);
	BOOST_CHECK(larger_subsequence_01[1] == Nucleotide::C);
	BOOST_CHECK(larger_subsequence_01[2] == Nucleotide::C);
	BOOST_CHECK(larger_subsequence_01[2] == Nucleotide::C);
	BOOST_CHECK(larger_subsequence_01.size == 4);
	
	
	// test the compare function
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_subsequence_01), IS_LARGER);// left sequence is larger than the 'smaller_' subsequence
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_subsequence_02), IS_LARGER);
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_subsequence_03), IS_LARGER);
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_subsequence_04), IS_LARGER);
	BOOST_REQUIRE_EQUAL(sequence.compare(smaller_subsequence_05), IS_LARGER);
	
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_subsequence_01), IS_SMALLER);// left sequence is larger than the 'larger_' subsequence
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_subsequence_02), IS_SMALLER);
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_subsequence_03), IS_SMALLER);
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_subsequence_04), IS_SMALLER);
	BOOST_REQUIRE_EQUAL(sequence.compare(larger_subsequence_05), IS_SMALLER);
}

BOOST_AUTO_TEST_SUITE_END()
