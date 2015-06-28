/**
 * @file test/Nucleotide_test.cpp
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



#define BOOST_TEST_MODULE Nucleotide



#include "Nucleotide.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether Nucleotides correctly separate
 *
 * @test
 *
 * @date 2014-04-15
 */
BOOST_AUTO_TEST_CASE(Test_equality)
{
	BOOST_CHECK_EQUAL(Nucleotide::A, Nucleotide::A);
	BOOST_CHECK_EQUAL(Nucleotide::C, Nucleotide::C);
	BOOST_CHECK_EQUAL(Nucleotide::G, Nucleotide::G);
	
	BOOST_CHECK_EQUAL(Nucleotide::U, Nucleotide::U);
	BOOST_CHECK_EQUAL(Nucleotide::T, Nucleotide::T);
	
	BOOST_CHECK_EQUAL(Nucleotide::U, Nucleotide::T);
	BOOST_CHECK_EQUAL(Nucleotide::T, Nucleotide::U);
	
	
	BOOST_CHECK(Nucleotide::A != Nucleotide::C);
	BOOST_CHECK(Nucleotide::A != Nucleotide::G);
	BOOST_CHECK(Nucleotide::A != Nucleotide::U);
	BOOST_CHECK(Nucleotide::A != Nucleotide::T);
	
	BOOST_CHECK(Nucleotide::C != Nucleotide::A);
	BOOST_CHECK(Nucleotide::C != Nucleotide::G);
	BOOST_CHECK(Nucleotide::C != Nucleotide::U);
	BOOST_CHECK(Nucleotide::C != Nucleotide::T);
	
	BOOST_CHECK(Nucleotide::G != Nucleotide::A);
	BOOST_CHECK(Nucleotide::G != Nucleotide::C);
	BOOST_CHECK(Nucleotide::G != Nucleotide::U);
	BOOST_CHECK(Nucleotide::G != Nucleotide::T);
	
	BOOST_CHECK(Nucleotide::U != Nucleotide::A);
	BOOST_CHECK(Nucleotide::T != Nucleotide::A);
	BOOST_CHECK(Nucleotide::U != Nucleotide::C);
	BOOST_CHECK(Nucleotide::T != Nucleotide::C);
	BOOST_CHECK(Nucleotide::U != Nucleotide::G);
	BOOST_CHECK(Nucleotide::T != Nucleotide::G);
}

/**
 * @brief Tests whether the size of a Nucleotide stuct is indeed 1 byte
 *
 * @test
 *
 * @date 2014-04-15
 */
BOOST_AUTO_TEST_CASE(Test_size)
{
	BOOST_CHECK_EQUAL(sizeof(Nucleotide::A), 1);
	BOOST_CHECK_EQUAL(sizeof(Nucleotide::C), 1);
	BOOST_CHECK_EQUAL(sizeof(Nucleotide::G), 1);
	
	BOOST_CHECK_EQUAL(sizeof(Nucleotide::U), 1);
	BOOST_CHECK_EQUAL(sizeof(Nucleotide::T), 1);
}

BOOST_AUTO_TEST_SUITE_END()
