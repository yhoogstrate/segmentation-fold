/**
 * @file test/Position_test.cpp
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



#define BOOST_TEST_MODULE Position



#include "main.hpp"
#include "Nucleotide.hpp"
#include "Position.hpp"
#include "Sequence.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether Nucleotides correctly separate
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	Sequence sequence_A = Sequence("A");
	Sequence sequence_C = Sequence("C");
	Sequence sequence_U = Sequence("U");
	Sequence sequence_G = Sequence("G");
	
	Position position_A = sequence_A.data.begin();
	Position position_C = sequence_C.data.begin();
	Position position_U = sequence_U.data.begin();
	Position position_G = sequence_G.data.begin();
	
	BOOST_CHECK_EQUAL((*position_A), Nucleotide::A);
	BOOST_CHECK_EQUAL((*position_C), Nucleotide::C);
	BOOST_CHECK_EQUAL((*position_U), Nucleotide::U);
	BOOST_CHECK_EQUAL((*position_U), Nucleotide::T);
	BOOST_CHECK_EQUAL((*position_G), Nucleotide::G);
}



/**
 * @brief Checks whether begin() last() and substractions with positions work
 *
 * @test
 *
 * @section DESCRIPTION
 * The following sequence "ACUG" will produce the following vector:
 * [Nucleotide::A][Nucleotide::C][Nucleotide::U][Nucleotide::G][ << end >> ]
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	Sequence sequence = Sequence("ACTG");
	
	Position position_begin = sequence.data.begin();
	Position position_last  = sequence.data.end() - 1;
	
	
	BOOST_CHECK_EQUAL(*position_begin, Nucleotide::A);
	BOOST_CHECK_EQUAL(*position_last,  Nucleotide::G);
	
	BOOST_CHECK_EQUAL(position_last - position_begin, 3);				// if normal integers would have been used, 3 - 0 = 3
}



/**
 * @brief check whether begin() last() and substractions with positions work and whether everything is fine after a sequence push_back
 *
 * @test
 *
 * @section DESCRIPTION
 * The following sequence "ACUG" will be the following vector:
 * [Nucleotide::A][Nucleotide::C][Nucleotide::U][Nucleotide::G][ << end >> ]
 *
 * Thus, a for loop works as follows:
 * for(Position p = s.begin(); p != s.end() ; p++) { ... }
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	Sequence sequence_poly_a = Sequence("AAAAAAAA");
	Sequence sequence_poly_c = Sequence("CccCCcCcc");
	Sequence sequence_poly_u = Sequence("UtTtu");
	Sequence sequence_poly_g = Sequence("G");
	Sequence sequence_empty  = Sequence("");
	
	Position i;
	for(i = sequence_poly_a.data.begin(); i != sequence_poly_a.data.end(); ++i)
	{
		BOOST_CHECK_EQUAL(*i, Nucleotide::A);
		
		BOOST_CHECK(*i != Nucleotide::C);
		BOOST_CHECK(*i != Nucleotide::T);
		BOOST_CHECK(*i != Nucleotide::U);
		BOOST_CHECK(*i != Nucleotide::G);
	}
	
	for(i = sequence_poly_c.data.begin(); i != sequence_poly_c.data.end(); ++i)
	{
		BOOST_CHECK(*i != Nucleotide::A);
		
		BOOST_CHECK_EQUAL(*i, Nucleotide::C);
		
		BOOST_CHECK(*i != Nucleotide::T);
		BOOST_CHECK(*i != Nucleotide::U);
		BOOST_CHECK(*i != Nucleotide::G);
	}
	
	for(i = sequence_poly_u.data.begin(); i != sequence_poly_u.data.end(); ++i)
	{
		BOOST_CHECK(*i != Nucleotide::A);
		BOOST_CHECK(*i != Nucleotide::C);
		
		BOOST_CHECK_EQUAL(*i, Nucleotide::T);
		BOOST_CHECK_EQUAL(*i, Nucleotide::U);
		
		BOOST_CHECK(*i != Nucleotide::G);
	}
	
	for(i = sequence_poly_g.data.begin(); i != sequence_poly_g.data.end(); ++i)
	{
		BOOST_CHECK(*i != Nucleotide::A);
		BOOST_CHECK(*i != Nucleotide::C);
		BOOST_CHECK(*i != Nucleotide::T);
		BOOST_CHECK(*i != Nucleotide::U);
		
		BOOST_CHECK_EQUAL(*i, Nucleotide::G);
	}
	
	for(i = sequence_empty.data.begin(); i != sequence_empty.data.end(); ++i)
	{
		// An empty sequence should never match any nucleotide - *i should be the << end >> element
		
		BOOST_CHECK(*i != Nucleotide::A);
		BOOST_CHECK(*i != Nucleotide::C);
		BOOST_CHECK(*i != Nucleotide::T);
		BOOST_CHECK(*i != Nucleotide::U);
		BOOST_CHECK(*i != Nucleotide::G);
	}
}



/**
 * @brief Checks whether end() and push_back() are in sync
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test4)
{
	Sequence sequence = Sequence();
	Position i;
	
	// Push back a Nucleotide::A
	sequence.push_back(Nucleotide::A);
	i = sequence.data.end() - 1;
	
	BOOST_CHECK_EQUAL(*i, Nucleotide::A);
	
	BOOST_CHECK(*i != Nucleotide::C);
	BOOST_CHECK(*i != Nucleotide::T);
	BOOST_CHECK(*i != Nucleotide::U);
	BOOST_CHECK(*i != Nucleotide::G);
	
	
	// Push back a Nucleotide::C
	sequence.push_back(Nucleotide::C);
	i = sequence.data.end() - 1;
	
	BOOST_CHECK(*i != Nucleotide::A);
	
	BOOST_CHECK_EQUAL(*i, Nucleotide::C);
	
	BOOST_CHECK(*i != Nucleotide::T);
	BOOST_CHECK(*i != Nucleotide::U);
	BOOST_CHECK(*i != Nucleotide::G);
	
	
	// Push back a Nucleotide::U
	sequence.push_back(Nucleotide::U);
	i = sequence.data.end() - 1;
	
	BOOST_CHECK(*i != Nucleotide::A);
	BOOST_CHECK(*i != Nucleotide::C);
	
	BOOST_CHECK_EQUAL(*i, Nucleotide::T);
	BOOST_CHECK_EQUAL(*i, Nucleotide::U);
	
	BOOST_CHECK(*i != Nucleotide::G);
	
	
	// Push back a Nucleotide::G
	sequence.push_back(Nucleotide::G);
	i = sequence.data.end() - 1;
	
	BOOST_CHECK(*i != Nucleotide::A);
	BOOST_CHECK(*i != Nucleotide::C);
	BOOST_CHECK(*i != Nucleotide::T);
	BOOST_CHECK(*i != Nucleotide::U);
	
	BOOST_CHECK_EQUAL(*i, Nucleotide::G);
}

BOOST_AUTO_TEST_SUITE_END()
