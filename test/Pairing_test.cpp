/**
 * @file test/Pairing_test.cpp
 *
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



#define BOOST_TEST_MODULE Pairing



#include "Nucleotide.hpp"
#include "Pairing.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether the is_canonical() function works
 *
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	Nucleotide a = Nucleotide::A;
	Nucleotide c = Nucleotide::C;
	Nucleotide u = Nucleotide::U;
	Nucleotide g = Nucleotide::G;
	
	Pairing aa = Pairing(a, a);
	Pairing ac = Pairing(a, c);
	Pairing au = Pairing(a, u);// Pair
	Pairing ag = Pairing(a, g);
	
	Pairing ca = Pairing(c, a);
	Pairing cc = Pairing(c, c);
	Pairing cu = Pairing(c, u);
	Pairing cg = Pairing(c, g);// Pair
	
	Pairing ua = Pairing(u, a);// Pair
	Pairing uc = Pairing(u, c);
	Pairing uu = Pairing(u, u);
	Pairing ug = Pairing(u, g);// Pair
	
	Pairing ga = Pairing(g, a);
	Pairing gc = Pairing(g, c);// Pair
	Pairing gu = Pairing(g, u);// Pair
	Pairing gg = Pairing(g, g);
	
	BOOST_CHECK(gc.is_canonical());
	BOOST_CHECK(cg.is_canonical());
	
	BOOST_CHECK(au.is_canonical());
	BOOST_CHECK(ua.is_canonical());
	
	// Wobbly pair
	BOOST_CHECK(gu.is_canonical());
	BOOST_CHECK(ug.is_canonical());
	
	
	BOOST_CHECK(aa.is_canonical() == false);
	BOOST_CHECK(ac.is_canonical() == false);
	BOOST_CHECK(ag.is_canonical() == false);
	
	BOOST_CHECK(ca.is_canonical() == false);
	BOOST_CHECK(cc.is_canonical() == false);
	BOOST_CHECK(cu.is_canonical() == false);
	
	BOOST_CHECK(uc.is_canonical() == false);
	BOOST_CHECK(uu.is_canonical() == false);
	
	BOOST_CHECK(ga.is_canonical() == false);
	BOOST_CHECK(gg.is_canonical() == false);
}

BOOST_AUTO_TEST_SUITE_END()
