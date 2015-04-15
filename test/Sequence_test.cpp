/**
 * @file test/Sequence_test.cpp
 *
 * @date 21-feb-2015
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



#define BOOST_TEST_MODULE Sequence



#include "Nucleotide.hpp"
#include "Sequence.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether single nucleotides are inserted correctly, and whether the size of the sequence fits.
 *
 * @date 11-mar-2014
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	Sequence s;
	
	s = Sequence("ACTG");
	BOOST_CHECK(s.size() == 4);
	
	Nucleotide n1 = Nucleotide::A;
	const char n2 = 'C';
	char n3 = 'u';
	
	s.push_back(n1);
	s.push_back(n2);
	s.push_back(n3);
	s.push_back('G');
	
	BOOST_CHECK_EQUAL(s.size() , 8);
	
	BOOST_CHECK_EQUAL(s[0] , Nucleotide::A);
	BOOST_CHECK_EQUAL(s[1] , Nucleotide::C);
	BOOST_CHECK_EQUAL(s[2] , Nucleotide::U);
	BOOST_CHECK_EQUAL(s[3] , Nucleotide::G);
	
	BOOST_CHECK_EQUAL(s[4] , Nucleotide::A);
	BOOST_CHECK_EQUAL(s[5] , Nucleotide::C);
	BOOST_CHECK_EQUAL(s[6] , Nucleotide::U);
	BOOST_CHECK_EQUAL(s[7] , Nucleotide::G);
}

/**
 * @brief Checks std::string to Sequence and Sequence to std::string conversion
 *
 * @date 12-mar-2014
 */
BOOST_AUTO_TEST_CASE(Test2)
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
	BOOST_CHECK(s2_str  .compare(s4.str()) == 0);
}

/**
 * @brief Tests behaviour of unexpected bases char-wise: should throw an exception
 *
 * @date 15-apr-2015
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	Sequence sequence;
	unsigned char c;
	
	for(c = 0; c < 255; c++)
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
 * @date 15-apr-2015
 */
BOOST_AUTO_TEST_CASE(Test4)
{
	Sequence sequence;
	std::string sequence_str;
	unsigned char c;
	
	for(c = 0; c < 255; c++)
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
 * @date 19-may-2014
 *
 * @section DESCRIPTION
 * Given that nucleotides are enumerated as follows:
 * A = 0
 * C = 1
 * G = 2
 * U = 3
 * A sequence of AA must be the 'smallest' of the two-letter sequences.
 */
BOOST_AUTO_TEST_CASE(Test5)
{
	Sequence s1 = Sequence("AA");
	Sequence s2 = Sequence("CC");
	Sequence s3 = Sequence("CG");
	Sequence s4a = Sequence("GU");
	Sequence s4b = Sequence("GT");
	Sequence s5 = Sequence("TU");
	
	
	BOOST_CHECK(s1 < s2);
	BOOST_CHECK(s2 < s3);
	BOOST_CHECK(s3 < s4a);
	BOOST_CHECK(s3 < s4b);
	BOOST_CHECK(s4a < s5);
	BOOST_CHECK(s4b < s5);
	
	BOOST_CHECK(s1 < s3);
	BOOST_CHECK(s2 < s4a);
	BOOST_CHECK(s2 < s4b);
	BOOST_CHECK(s3 < s5);
	
	BOOST_CHECK(s1 < s4a);
	BOOST_CHECK(s1 < s4b);
	BOOST_CHECK(s2 < s5);
	
	BOOST_CHECK(s1 < s5);
	
	
	BOOST_CHECK((s2 < s1) == false);
	BOOST_CHECK((s3 < s2) == false);
	BOOST_CHECK((s4a < s3) == false);
	BOOST_CHECK((s4b < s3) == false);
	BOOST_CHECK((s5 < s4a) == false);
	BOOST_CHECK((s5 < s4b) == false);
	
	BOOST_CHECK((s3 < s1) == false);
	BOOST_CHECK((s4a < s2) == false);
	BOOST_CHECK((s4b < s2) == false);
	BOOST_CHECK((s5 < s3) == false);
	
	BOOST_CHECK((s4a < s1) == false);
	BOOST_CHECK((s4b < s1) == false);
	BOOST_CHECK((s5 < s2) == false);
	
	BOOST_CHECK((s5 < s1) == false);
}

/**
 * @brief Tests the '==' operator
 *
 * @date 19-may-2014
 */
BOOST_AUTO_TEST_CASE(Test6)
{
	Sequence s1 = Sequence("AAGAA");
	Sequence s2 = Sequence("aAGAA");
	Sequence s3 = Sequence("AAGA");
	Sequence s4 = Sequence("AAGAAA");
	Sequence s5 = Sequence("TTGTT");
	Sequence s6 = Sequence("TUGuu");
	
	
	BOOST_CHECK(s1 == s2);
	BOOST_CHECK((s2 == s3) == false);
	BOOST_CHECK((s3 == s4) == false);
	BOOST_CHECK((s4 == s5) == false);
	BOOST_CHECK(s5 == s6);
	
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
 *
 * @date 19-feb-2015
 */
BOOST_AUTO_TEST_CASE(Test7)
{
	Sequence sequence = Sequence("ACTG");
	
	BOOST_CHECK_EQUAL(sequence.size() , 4);
	
	
	Sequence subsequence_01 = sequence.subseq(0, 0);// start at 0, stop at 0 -> A
	Sequence subsequence_02 = sequence.subseq(0, 1);// start at 0, stop at 1 -> AC
	Sequence subsequence_03 = sequence.subseq(0, 2);// start at 0, stop at 2 -> ACT
	Sequence subsequence_04 = sequence.subseq(0, 3);// start at 0, stop at 3 -> ACTG
	Sequence subsequence_05 = sequence.subseq(1, 3);// start at 1, stop at 3 ->  CTG
	Sequence subsequence_06 = sequence.subseq(2, 3);// start at 2, stop at 3 ->   TG
	Sequence subsequence_07 = sequence.subseq(3, 3);// start at 3, stop at 3 ->    G
	
	BOOST_CHECK_EQUAL(subsequence_01.size() , 1);
	BOOST_CHECK_EQUAL(subsequence_02.size() , 2);
	BOOST_CHECK_EQUAL(subsequence_03.size() , 3);
	BOOST_CHECK_EQUAL(subsequence_04.size() , 4);
	BOOST_CHECK_EQUAL(subsequence_05.size() , 3);
	BOOST_CHECK_EQUAL(subsequence_06.size() , 2);
	BOOST_CHECK_EQUAL(subsequence_07.size() , 1);
	
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
BOOST_AUTO_TEST_SUITE_END()
