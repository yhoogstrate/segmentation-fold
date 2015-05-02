/**
 * @file test/PairingPlus_test.cpp
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



#define BOOST_TEST_MODULE PairingPlus



#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "PairingType.hpp"
#include "PairingPlus.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether is_canonical() and initialization works
 * 
 * @test
 * 
 * @date 2015-05-02
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	Sequence sequence = Sequence("ACUG");
	Position a = sequence.data.begin();
	Position c = sequence.data.begin() + 1;
	Position u = sequence.data.begin() + 2;
	Position g = sequence.data.begin() + 3;
	
	PairingPlus aa = PairingPlus(a, a);
	PairingPlus ac = PairingPlus(a, c);
	PairingPlus au = PairingPlus(a, u);// True canonical pair
	PairingPlus ag = PairingPlus(a, g);
	
	PairingPlus ca = PairingPlus(c, a);
	PairingPlus cc = PairingPlus(c, c);
	PairingPlus cu = PairingPlus(c, u);
	PairingPlus cg = PairingPlus(c, g);// True canonical pair
	
	PairingPlus ua = PairingPlus(u, a);// True canonical pair
	PairingPlus uc = PairingPlus(u, c);
	PairingPlus uu = PairingPlus(u, u);
	PairingPlus ug = PairingPlus(u, g);// True canonical pair
	
	PairingPlus ga = PairingPlus(g, a);
	PairingPlus gc = PairingPlus(g, c);// True canonical pair
	PairingPlus gu = PairingPlus(g, u);// True canonical pair
	PairingPlus gg = PairingPlus(g, g);
	
	BOOST_CHECK(gc.is_canonical());
	BOOST_CHECK(cg.is_canonical());
	
	BOOST_CHECK(au.is_canonical());
	BOOST_CHECK(ua.is_canonical());
	
	// Wobbly pairs
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
	
	
	
	BOOST_CHECK(aa.type == PairingType::None);
	BOOST_CHECK(ac.type == PairingType::None);
	BOOST_CHECK(au.type == PairingType::AU);// True canonical pair
	BOOST_CHECK(ag.type == PairingType::None);
	
	BOOST_CHECK(ca.type == PairingType::None);
	BOOST_CHECK(cc.type == PairingType::None);
	BOOST_CHECK(cu.type == PairingType::None);
	BOOST_CHECK(cg.type == PairingType::CG);// True canonical pair
	
	BOOST_CHECK(ua.type == PairingType::UA);// True canonical pair
	BOOST_CHECK(uc.type == PairingType::None);
	BOOST_CHECK(uu.type == PairingType::None);
	BOOST_CHECK(ug.type == PairingType::UG);// True canonical pair
	
	BOOST_CHECK(ga.type == PairingType::None);
	BOOST_CHECK(gc.type == PairingType::GC);// True canonical pair
	BOOST_CHECK(gu.type == PairingType::GU);// True canonical pair
	BOOST_CHECK(gg.type == PairingType::None);
}

BOOST_AUTO_TEST_SUITE_END()
