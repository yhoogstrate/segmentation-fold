/**
 * @file test/Direction_test.cpp
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



#define BOOST_TEST_MODULE Direction



#include "Direction.hpp"
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether 5' and 3' are not identical
 *
 * @test
 *
 * @date 2015-04-23
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	BOOST_CHECK(Direction::FivePrime  == Direction::FivePrime);
	BOOST_CHECK(Direction::ThreePrime == Direction::ThreePrime);
	BOOST_CHECK(Direction::FivePrime  != Direction::ThreePrime);
	
	
	Direction d5_1 = Direction::FivePrime;
	Direction d5_2 = Direction::FivePrime;
	
	Direction d3_1 = Direction::ThreePrime;
	Direction d3_2 = Direction::ThreePrime;
	
	BOOST_CHECK(d5_1 == d5_2);
	BOOST_CHECK(d3_1 == d3_2);
	
	BOOST_CHECK(d5_1 != d3_1);
	BOOST_CHECK(d5_1 != d3_2);
	
	BOOST_CHECK(d5_2 != d3_1);
	BOOST_CHECK(d5_2 != d3_2);
}



/**
 * @brief Tests whether the size of a direction is 1 byte
 *
 * @test
 *
 * @date 2015-04-22
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	Direction d1 = Direction::FivePrime;
	Direction d2 = Direction::ThreePrime;
	
	BOOST_CHECK_EQUAL(sizeof(d1), 1);
	BOOST_CHECK_EQUAL(sizeof(d2), 1);
}

BOOST_AUTO_TEST_SUITE_END()
