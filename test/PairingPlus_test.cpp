/**
 * @file test/PairingPlus_test.cpp
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



#define BOOST_TEST_MODULE PairingPlus



#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "PairingType.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether is_canonical() and initialization works
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	Sequence sequence = Sequence("ACUGACUG");
	Position a1 = sequence.data.begin() + 0;
	Position c1 = sequence.data.begin() + 1;
	Position u1 = sequence.data.begin() + 2;
	Position g1 = sequence.data.begin() + 3;
	Position a2 = sequence.data.begin() + 4;
	Position c2 = sequence.data.begin() + 5;
	Position u2 = sequence.data.begin() + 6;
	Position g2 = sequence.data.begin() + 7;
	
	PairingPlus aa = PairingPlus(a1, a2);
	PairingPlus ac = PairingPlus(a1, c2);
	PairingPlus au = PairingPlus(a1, u2);// True canonical pair
	PairingPlus ag = PairingPlus(a1, g2);
	
	PairingPlus ca = PairingPlus(c1, a2);
	PairingPlus cc = PairingPlus(c1, c2);
	PairingPlus cu = PairingPlus(c1, u2);
	PairingPlus cg = PairingPlus(c1, g2);// True canonical pair
	
	PairingPlus ua = PairingPlus(u1, a2);// True canonical pair
	PairingPlus uc = PairingPlus(u1, c2);
	PairingPlus uu = PairingPlus(u1, u2);
	PairingPlus ug = PairingPlus(u1, g2);// True canonical pair
	
	PairingPlus ga = PairingPlus(g1, a2);
	PairingPlus gc = PairingPlus(g1, c2);// True canonical pair
	PairingPlus gu = PairingPlus(g1, u2);// True canonical pair
	PairingPlus gg = PairingPlus(g1, g2);
	
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
	
	
	
	BOOST_CHECK_EQUAL(aa.type , PairingType::None);
	BOOST_CHECK_EQUAL(ac.type , PairingType::None);
	BOOST_CHECK_EQUAL(au.type , PairingType::AU);// True canonical pair
	BOOST_CHECK_EQUAL(ag.type , PairingType::None);
	
	BOOST_CHECK_EQUAL(ca.type , PairingType::None);
	BOOST_CHECK_EQUAL(cc.type , PairingType::None);
	BOOST_CHECK_EQUAL(cu.type , PairingType::None);
	BOOST_CHECK_EQUAL(cg.type , PairingType::CG);// True canonical pair
	
	BOOST_CHECK_EQUAL(ua.type , PairingType::UA);// True canonical pair
	BOOST_CHECK_EQUAL(uc.type , PairingType::None);
	BOOST_CHECK_EQUAL(uu.type , PairingType::None);
	BOOST_CHECK_EQUAL(ug.type , PairingType::UG);// True canonical pair
	
	BOOST_CHECK_EQUAL(ga.type , PairingType::None);
	BOOST_CHECK_EQUAL(gc.type , PairingType::GC);// True canonical pair
	BOOST_CHECK_EQUAL(gu.type , PairingType::GU);// True canonical pair
	BOOST_CHECK_EQUAL(gg.type , PairingType::None);
	
	
	// Size is size in between the pairs:
	// CA
	// | A
	// GA
	// Has size 3
	BOOST_CHECK_EQUAL(aa.size , 3);
	BOOST_CHECK_EQUAL(ac.size , 4);
	BOOST_CHECK_EQUAL(au.size , 5);
	BOOST_CHECK_EQUAL(ag.size , 6);
	
	BOOST_CHECK_EQUAL(ca.size , 2);
	BOOST_CHECK_EQUAL(cc.size , 3);
	BOOST_CHECK_EQUAL(cu.size , 4);
	BOOST_CHECK_EQUAL(cg.size , 5);
	
	BOOST_CHECK_EQUAL(ua.size , 1);
	BOOST_CHECK_EQUAL(uc.size , 2);
	BOOST_CHECK_EQUAL(uu.size , 3);
	BOOST_CHECK_EQUAL(ug.size , 4);
	
	BOOST_CHECK_EQUAL(ga.size , 0);
	BOOST_CHECK_EQUAL(gc.size , 1);
	BOOST_CHECK_EQUAL(gu.size , 2);
	BOOST_CHECK_EQUAL(gg.size , 3);
}



/**
 * @brief Tests whether is_canonical() and initialization using Sequence works
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	Sequence sequence = Sequence("ACUGACUG");
	char a1 = 0;
	char c1 = 1;
	char u1 = 2;
	char g1 = 3;
	char a2 = 4;
	char c2 = 5;
	char u2 = 6;
	char g2 = 7;
	
	Pair paa = Pair(a1, a2);
	Pair pac = Pair(a1, c2);
	Pair pau = Pair(a1, u2);
	Pair pag = Pair(a1, g2);
	
	Pair pca = Pair(c1, a2);
	Pair pcc = Pair(c1, c2);
	Pair pcu = Pair(c1, u2);
	Pair pcg = Pair(c1, g2);
	
	Pair pua = Pair(u1, a2);
	Pair puc = Pair(u1, c2);
	Pair puu = Pair(u1, u2);
	Pair pug = Pair(u1, g2);
	
	Pair pga = Pair(g1, a2);
	Pair pgc = Pair(g1, c2);
	Pair pgu = Pair(g1, u2);
	Pair pgg = Pair(g1, g2);
	
	
	
	PairingPlus aa = PairingPlus(sequence, paa);
	PairingPlus ac = PairingPlus(sequence, pac);
	PairingPlus au = PairingPlus(sequence, pau);
	PairingPlus ag = PairingPlus(sequence, pag);
	
	PairingPlus ca = PairingPlus(sequence, pca);
	PairingPlus cc = PairingPlus(sequence, pcc);
	PairingPlus cu = PairingPlus(sequence, pcu);
	PairingPlus cg = PairingPlus(sequence, pcg);
	
	PairingPlus ua = PairingPlus(sequence, pua);
	PairingPlus uc = PairingPlus(sequence, puc);
	PairingPlus uu = PairingPlus(sequence, puu);
	PairingPlus ug = PairingPlus(sequence, pug);
	
	PairingPlus ga = PairingPlus(sequence, pga);
	PairingPlus gc = PairingPlus(sequence, pgc);
	PairingPlus gu = PairingPlus(sequence, pgu);
	PairingPlus gg = PairingPlus(sequence, pgg);
	
	
	
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
	
	
	
	BOOST_CHECK_EQUAL(aa.type , PairingType::None);
	BOOST_CHECK_EQUAL(ac.type , PairingType::None);
	BOOST_CHECK_EQUAL(au.type , PairingType::AU);// True canonical pair
	BOOST_CHECK_EQUAL(ag.type , PairingType::None);
	
	BOOST_CHECK_EQUAL(ca.type , PairingType::None);
	BOOST_CHECK_EQUAL(cc.type , PairingType::None);
	BOOST_CHECK_EQUAL(cu.type , PairingType::None);
	BOOST_CHECK_EQUAL(cg.type , PairingType::CG);// True canonical pair
	
	BOOST_CHECK_EQUAL(ua.type , PairingType::UA);// True canonical pair
	BOOST_CHECK_EQUAL(uc.type , PairingType::None);
	BOOST_CHECK_EQUAL(uu.type , PairingType::None);
	BOOST_CHECK_EQUAL(ug.type , PairingType::UG);// True canonical pair
	
	BOOST_CHECK_EQUAL(ga.type , PairingType::None);
	BOOST_CHECK_EQUAL(gc.type , PairingType::GC);// True canonical pair
	BOOST_CHECK_EQUAL(gu.type , PairingType::GU);// True canonical pair
	BOOST_CHECK_EQUAL(gg.type , PairingType::None);
	
	
	// Size is size in between the pairs:
	// CA
	// | A
	// GA
	// Has size 3
	BOOST_CHECK_EQUAL(aa.size , 3);
	BOOST_CHECK_EQUAL(ac.size , 4);
	BOOST_CHECK_EQUAL(au.size , 5);
	BOOST_CHECK_EQUAL(ag.size , 6);
	
	BOOST_CHECK_EQUAL(ca.size , 2);
	BOOST_CHECK_EQUAL(cc.size , 3);
	BOOST_CHECK_EQUAL(cu.size , 4);
	BOOST_CHECK_EQUAL(cg.size , 5);
	
	BOOST_CHECK_EQUAL(ua.size , 1);
	BOOST_CHECK_EQUAL(uc.size , 2);
	BOOST_CHECK_EQUAL(uu.size , 3);
	BOOST_CHECK_EQUAL(ug.size , 4);
	
	BOOST_CHECK_EQUAL(ga.size , 0);
	BOOST_CHECK_EQUAL(gc.size , 1);
	BOOST_CHECK_EQUAL(gu.size , 2);
	BOOST_CHECK_EQUAL(gg.size , 3);
}

BOOST_AUTO_TEST_SUITE_END()
