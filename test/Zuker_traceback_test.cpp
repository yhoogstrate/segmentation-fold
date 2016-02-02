/**
 * @file test/Zuker_traceback_test.cpp
 *
 * @date 2016-01-22
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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
 * </PRE>
 */



#define BOOST_TEST_MODULE main


#include "main.hpp"

#include "Pair.hpp"
#include "Region.hpp"
#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"
#include "SegmentLoop.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"
#include "SegmentLoopTree.hpp"
#include "ReadSegments.hpp"

#include "ScoringMatrix.hpp"
#include "Settings.hpp"
#include "DotBracket.hpp"
#include "ReadData.hpp"

#include "Zuker.hpp"
#include "ZukerTraceback.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)


bool IdenticalFloats(double a, double b)
{
	return fabs(a - b) < 0.0001;
}



/**
 * @date 2016-01-22
 */
void test_vij(Zuker &zuker)
{
	Pair p = Pair();
	float e;
	
	p = {0, 0};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(0,0) = " + std::to_string(e));
	p = {0, 1};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(0,1) = " + std::to_string(e));
	p = {0, 2};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(0,2) = " + std::to_string(e));
	p = {0, 3};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(0,3) = " + std::to_string(e));
	p = {0, 4};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,4) = " + std::to_string(e));
	p = {0, 5};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,5) = " + std::to_string(e));
	p = {0, 6};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,6) = " + std::to_string(e));
	p = {0, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,7) = " + std::to_string(e));
	p = {0, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,8) = " + std::to_string(e));
	p = {0, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,9) = " + std::to_string(e));
	p = {0, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,10) = " + std::to_string(e));
	p = {0, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,11) = " + std::to_string(e));
	p = {0, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(0,12) = " + std::to_string(e));
	p = {0, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(0,13) = " + std::to_string(e));
	p = {0, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(0,14) = " + std::to_string(e));
	p = {0, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(0,15) = " + std::to_string(e));
	p = {0, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.90f), "V(0,16) = " + std::to_string(e));
	p = {0, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.90f), "V(0,17) = " + std::to_string(e));
	
	p = {1, 1};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(1,1) = " + std::to_string(e));
	p = {1, 2};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(1,2) = " + std::to_string(e));
	p = {1, 3};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(1,3) = " + std::to_string(e));
	p = {1, 4};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(1,4) = " + std::to_string(e));
	p = {1, 5};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,5) = " + std::to_string(e));
	p = {1, 6};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,6) = " + std::to_string(e));
	p = {1, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,7) = " + std::to_string(e));
	p = {1, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,8) = " + std::to_string(e));
	p = {1, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,9) = " + std::to_string(e));
	p = {1, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,10) = " + std::to_string(e));
	p = {1, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,11) = " + std::to_string(e));
	p = {1, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(1,12) = " + std::to_string(e));
	p = {1, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(1,13) = " + std::to_string(e));
	p = {1, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(1,14) = " + std::to_string(e));
	p = {1, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(1,15) = " + std::to_string(e));
	p = {1, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(1,16) = " + std::to_string(e));
	p = {1, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.90f), "V(1,17) = " + std::to_string(e));
	
	p = {2, 2};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(2,2) = " + std::to_string(e));
	p = {2, 3};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(2,3) = " + std::to_string(e));
	p = {2, 4};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(2,4) = " + std::to_string(e));
	p = {2, 5};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(2,5) = " + std::to_string(e));
	p = {2, 6};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,6) = " + std::to_string(e));
	p = {2, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,7) = " + std::to_string(e));
	p = {2, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,8) = " + std::to_string(e));
	p = {2, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,9) = " + std::to_string(e));
	p = {2, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,10) = " + std::to_string(e));
	p = {2, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,11) = " + std::to_string(e));
	p = {2, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,12) = " + std::to_string(e));
	p = {2, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(2,13) = " + std::to_string(e));
	p = {2, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(2,14) = " + std::to_string(e));
	p = {2, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(2,15) = " + std::to_string(e));
	p = {2, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(2,16) = " + std::to_string(e));
	p = {2, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(2,17) = " + std::to_string(e));
	
	p = {3, 3};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(3,3) = " + std::to_string(e));
	p = {3, 4};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(3,4) = " + std::to_string(e));
	p = {3, 5};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(3,5) = " + std::to_string(e));
	p = {3, 6};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(3,6) = " + std::to_string(e));
	p = {3, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,7) = " + std::to_string(e));
	p = {3, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,8) = " + std::to_string(e));
	p = {3, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,9) = " + std::to_string(e));
	p = {3, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,10) = " + std::to_string(e));
	p = {3, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,11) = " + std::to_string(e));
	p = {3, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,12) = " + std::to_string(e));
	p = {3, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,13) = " + std::to_string(e));
	p = {3, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,14) = " + std::to_string(e));
	p = {3, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,15) = " + std::to_string(e));
	p = {3, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,16) = " + std::to_string(e));
	p = {3, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(3,17) = " + std::to_string(e));
	
	p = {4, 4};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(4,4) = " + std::to_string(e));
	p = {4, 5};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(4,5) = " + std::to_string(e));
	p = {4, 6};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(4,6) = " + std::to_string(e));
	p = {4, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(4,7) = " + std::to_string(e));
	p = {4, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,8) = " + std::to_string(e));
	p = {4, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,9) = " + std::to_string(e));
	p = {4, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,10) = " + std::to_string(e));
	p = {4, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,11) = " + std::to_string(e));
	p = {4, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,12) = " + std::to_string(e));
	p = {4, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,13) = " + std::to_string(e));
	p = {4, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,14) = " + std::to_string(e));
	p = {4, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,15) = " + std::to_string(e));
	p = {4, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,16) = " + std::to_string(e));
	p = {4, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(4,17) = " + std::to_string(e));
	
	p = {5, 5};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(5,5) = " + std::to_string(e));
	p = {5, 6};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(5,6) = " + std::to_string(e));
	p = {5, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(5,7) = " + std::to_string(e));
	p = {5, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(5,8) = " + std::to_string(e));
	p = {5, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,9) = " + std::to_string(e));
	p = {5, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,10) = " + std::to_string(e));
	p = {5, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,11) = " + std::to_string(e));
	p = {5, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,12) = " + std::to_string(e));
	p = {5, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,13) = " + std::to_string(e));
	p = {5, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,14) = " + std::to_string(e));
	p = {5, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,15) = " + std::to_string(e));
	p = {5, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,16) = " + std::to_string(e));
	p = {5, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(5,17) = " + std::to_string(e));
	
	p = {6, 6};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(6,6) = " + std::to_string(e));
	p = {6, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(6,7) = " + std::to_string(e));
	p = {6, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(6,8) = " + std::to_string(e));
	p = {6, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(6,9) = " + std::to_string(e));
	p = {6, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(6,10) = " + std::to_string(e));
	p = {6, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(6,11) = " + std::to_string(e));
	p = {6, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(6,12) = " + std::to_string(e));
	p = {6, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(6,13) = " + std::to_string(e));
	p = {6, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(6,14) = " + std::to_string(e));
	p = {6, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(6,15) = " + std::to_string(e));
	p = {6, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(6,16) = " + std::to_string(e));
	p = {6, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "V(6,17) = " + std::to_string(e));
	
	p = {7, 7};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(7,7) = " + std::to_string(e));
	p = {7, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(7,8) = " + std::to_string(e));
	p = {7, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(7,9) = " + std::to_string(e));
	p = {7, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(7,10) = " + std::to_string(e));
	p = {7, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(7,11) = " + std::to_string(e));
	p = {7, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(7,12) = " + std::to_string(e));
	p = {7, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(7,13) = " + std::to_string(e));
	p = {7, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(7,14) = " + std::to_string(e));
	p = {7, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(7,15) = " + std::to_string(e));
	p = {7, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(7,16) = " + std::to_string(e));
	p = {7, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "V(7,17) = " + std::to_string(e));
	
	p = {8, 8};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(8,8) = " + std::to_string(e));
	p = {8, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(8,9) = " + std::to_string(e));
	p = {8, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(8,10) = " + std::to_string(e));
	p = {8, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(8,11) = " + std::to_string(e));
	p = {8, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 4.60f), "V(8,12) = " + std::to_string(e));
	p = {8, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(8,13) = " + std::to_string(e));
	p = {8, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(8,14) = " + std::to_string(e));
	p = {8, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(8,15) = " + std::to_string(e));
	p = {8, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(8,16) = " + std::to_string(e));
	p = {8, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(8,17) = " + std::to_string(e));
	
	p = {9, 9};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(9,9) = " + std::to_string(e));
	p = {9, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(9,10) = " + std::to_string(e));
	p = {9, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(9,11) = " + std::to_string(e));
	p = {9, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(9,12) = " + std::to_string(e));
	p = {9, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(9,13) = " + std::to_string(e));
	p = {9, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(9,14) = " + std::to_string(e));
	p = {9, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(9,15) = " + std::to_string(e));
	p = {9, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(9,16) = " + std::to_string(e));
	p = {9, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(9,17) = " + std::to_string(e));
	
	p = {10, 10};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(10,10) = " + std::to_string(e));
	p = {10, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(10,11) = " + std::to_string(e));
	p = {10, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(10,12) = " + std::to_string(e));
	p = {10, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(10,13) = " + std::to_string(e));
	p = {10, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(10,14) = " + std::to_string(e));
	p = {10, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(10,15) = " + std::to_string(e));
	p = {10, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(10,16) = " + std::to_string(e));
	p = {10, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(10,17) = " + std::to_string(e));
	
	p = {11, 11};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(11,11) = " + std::to_string(e));
	p = {11, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(11,12) = " + std::to_string(e));
	p = {11, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(11,13) = " + std::to_string(e));
	p = {11, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(11,14) = " + std::to_string(e));
	p = {11, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(11,15) = " + std::to_string(e));
	p = {11, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(11,16) = " + std::to_string(e));
	p = {11, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(11,17) = " + std::to_string(e));
	
	p = {12, 12};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(12,12) = " + std::to_string(e));
	p = {12, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(12,13) = " + std::to_string(e));
	p = {12, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(12,14) = " + std::to_string(e));
	p = {12, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(12,15) = " + std::to_string(e));
	p = {12, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(12,16) = " + std::to_string(e));
	p = {12, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(12,17) = " + std::to_string(e));
	
	p = {13, 13};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(13,13) = " + std::to_string(e));
	p = {13, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(13,14) = " + std::to_string(e));
	p = {13, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(13,15) = " + std::to_string(e));
	p = {13, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(13,16) = " + std::to_string(e));
	p = {13, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "V(13,17) = " + std::to_string(e));
	
	p = {14, 14};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(14,14) = " + std::to_string(e));
	p = {14, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(14,15) = " + std::to_string(e));
	p = {14, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(14,16) = " + std::to_string(e));
	p = {14, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(14,17) = " + std::to_string(e));
	
	p = {15, 15};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(15,15) = " + std::to_string(e));
	p = {15, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(15,16) = " + std::to_string(e));
	p = {15, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(15,17) = " + std::to_string(e));
	
	p = {16, 16};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(16,16) = " + std::to_string(e));
	p = {16, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(16,17) = " + std::to_string(e));
	
	p = {17, 17};
	e = zuker.vij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, N_INFINITY), "V(17,17) = " + std::to_string(e));
}



/**
 * @date 2016-01-22
 */
void test_wij(Zuker &zuker)
{
	Pair p = Pair();
	float e;
	
	p = {0, 0};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,0) = " + std::to_string(e));
	p = {0, 1};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,1) = " + std::to_string(e));
	p = {0, 2};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,2) = " + std::to_string(e));
	p = {0, 3};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,3) = " + std::to_string(e));
	p = {0, 4};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,4) = " + std::to_string(e));
	p = {0, 5};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,5) = " + std::to_string(e));
	p = {0, 6};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,6) = " + std::to_string(e));
	p = {0, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,7) = " + std::to_string(e));
	p = {0, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,8) = " + std::to_string(e));
	p = {0, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,9) = " + std::to_string(e));
	p = {0, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,10) = " + std::to_string(e));
	p = {0, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,11) = " + std::to_string(e));
	p = {0, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(0,12) = " + std::to_string(e));
	p = {0, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(0,13) = " + std::to_string(e));
	p = {0, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(0,14) = " + std::to_string(e));
	p = {0, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(0,15) = " + std::to_string(e));
	p = {0, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.90f), "W(0,16) = " + std::to_string(e));
	p = {0, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.90f), "W(0,17) = " + std::to_string(e));
	
	p = {1, 1};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,1) = " + std::to_string(e));
	p = {1, 2};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,2) = " + std::to_string(e));
	p = {1, 3};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,3) = " + std::to_string(e));
	p = {1, 4};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,4) = " + std::to_string(e));
	p = {1, 5};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,5) = " + std::to_string(e));
	p = {1, 6};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,6) = " + std::to_string(e));
	p = {1, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,7) = " + std::to_string(e));
	p = {1, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,8) = " + std::to_string(e));
	p = {1, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,9) = " + std::to_string(e));
	p = {1, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,10) = " + std::to_string(e));
	p = {1, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,11) = " + std::to_string(e));
	p = {1, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(1,12) = " + std::to_string(e));
	p = {1, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(1,13) = " + std::to_string(e));
	p = {1, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(1,14) = " + std::to_string(e));
	p = {1, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(1,15) = " + std::to_string(e));
	p = {1, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(1,16) = " + std::to_string(e));
	p = {1, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.90f), "W(1,17) = " + std::to_string(e));
	
	p = {2, 2};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,2) = " + std::to_string(e));
	p = {2, 3};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,3) = " + std::to_string(e));
	p = {2, 4};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,4) = " + std::to_string(e));
	p = {2, 5};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,5) = " + std::to_string(e));
	p = {2, 6};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,6) = " + std::to_string(e));
	p = {2, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,7) = " + std::to_string(e));
	p = {2, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,8) = " + std::to_string(e));
	p = {2, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,9) = " + std::to_string(e));
	p = {2, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,10) = " + std::to_string(e));
	p = {2, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,11) = " + std::to_string(e));
	p = {2, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(2,12) = " + std::to_string(e));
	p = {2, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(2,13) = " + std::to_string(e));
	p = {2, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(2,14) = " + std::to_string(e));
	p = {2, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(2,15) = " + std::to_string(e));
	p = {2, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(2,16) = " + std::to_string(e));
	p = {2, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(2,17) = " + std::to_string(e));
	
	p = {3, 3};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,3) = " + std::to_string(e));
	p = {3, 4};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,4) = " + std::to_string(e));
	p = {3, 5};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,5) = " + std::to_string(e));
	p = {3, 6};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,6) = " + std::to_string(e));
	p = {3, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,7) = " + std::to_string(e));
	p = {3, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,8) = " + std::to_string(e));
	p = {3, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,9) = " + std::to_string(e));
	p = {3, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,10) = " + std::to_string(e));
	p = {3, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,11) = " + std::to_string(e));
	p = {3, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(3,12) = " + std::to_string(e));
	p = {3, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(3,13) = " + std::to_string(e));
	p = {3, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(3,14) = " + std::to_string(e));
	p = {3, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(3,15) = " + std::to_string(e));
	p = {3, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(3,16) = " + std::to_string(e));
	p = {3, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(3,17) = " + std::to_string(e));
	
	p = {4, 4};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,4) = " + std::to_string(e));
	p = {4, 5};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,5) = " + std::to_string(e));
	p = {4, 6};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,6) = " + std::to_string(e));
	p = {4, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,7) = " + std::to_string(e));
	p = {4, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,8) = " + std::to_string(e));
	p = {4, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,9) = " + std::to_string(e));
	p = {4, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,10) = " + std::to_string(e));
	p = {4, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,11) = " + std::to_string(e));
	p = {4, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(4,12) = " + std::to_string(e));
	p = {4, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(4,13) = " + std::to_string(e));
	p = {4, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(4,14) = " + std::to_string(e));
	p = {4, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(4,15) = " + std::to_string(e));
	p = {4, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(4,16) = " + std::to_string(e));
	p = {4, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(4,17) = " + std::to_string(e));
	
	p = {5, 5};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,5) = " + std::to_string(e));
	p = {5, 6};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,6) = " + std::to_string(e));
	p = {5, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,7) = " + std::to_string(e));
	p = {5, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,8) = " + std::to_string(e));
	p = {5, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,9) = " + std::to_string(e));
	p = {5, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,10) = " + std::to_string(e));
	p = {5, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,11) = " + std::to_string(e));
	p = {5, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(5,12) = " + std::to_string(e));
	p = {5, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(5,13) = " + std::to_string(e));
	p = {5, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(5,14) = " + std::to_string(e));
	p = {5, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(5,15) = " + std::to_string(e));
	p = {5, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(5,16) = " + std::to_string(e));
	p = {5, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(5,17) = " + std::to_string(e));
	
	p = {6, 6};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(6,6) = " + std::to_string(e));
	p = {6, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(6,7) = " + std::to_string(e));
	p = {6, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(6,8) = " + std::to_string(e));
	p = {6, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(6,9) = " + std::to_string(e));
	p = {6, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(6,10) = " + std::to_string(e));
	p = {6, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(6,11) = " + std::to_string(e));
	p = {6, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(6,12) = " + std::to_string(e));
	p = {6, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(6,13) = " + std::to_string(e));
	p = {6, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(6,14) = " + std::to_string(e));
	p = {6, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(6,15) = " + std::to_string(e));
	p = {6, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(6,16) = " + std::to_string(e));
	p = {6, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.60f), "W(6,17) = " + std::to_string(e));
	
	p = {7, 7};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(7,7) = " + std::to_string(e));
	p = {7, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(7,8) = " + std::to_string(e));
	p = {7, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(7,9) = " + std::to_string(e));
	p = {7, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(7,10) = " + std::to_string(e));
	p = {7, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(7,11) = " + std::to_string(e));
	p = {7, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(7,12) = " + std::to_string(e));
	p = {7, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(7,13) = " + std::to_string(e));
	p = {7, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(7,14) = " + std::to_string(e));
	p = {7, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(7,15) = " + std::to_string(e));
	p = {7, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(7,16) = " + std::to_string(e));
	p = {7, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.30f), "W(7,17) = " + std::to_string(e));
	
	p = {8, 8};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,8) = " + std::to_string(e));
	p = {8, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,9) = " + std::to_string(e));
	p = {8, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,10) = " + std::to_string(e));
	p = {8, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,11) = " + std::to_string(e));
	p = {8, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,12) = " + std::to_string(e));
	p = {8, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,13) = " + std::to_string(e));
	p = {8, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,14) = " + std::to_string(e));
	p = {8, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,15) = " + std::to_string(e));
	p = {8, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,16) = " + std::to_string(e));
	p = {8, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(8,17) = " + std::to_string(e));
	
	p = {9, 9};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,9) = " + std::to_string(e));
	p = {9, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,10) = " + std::to_string(e));
	p = {9, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,11) = " + std::to_string(e));
	p = {9, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,12) = " + std::to_string(e));
	p = {9, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,13) = " + std::to_string(e));
	p = {9, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,14) = " + std::to_string(e));
	p = {9, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,15) = " + std::to_string(e));
	p = {9, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,16) = " + std::to_string(e));
	p = {9, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(9,17) = " + std::to_string(e));
	
	p = {10, 10};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,10) = " + std::to_string(e));
	p = {10, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,11) = " + std::to_string(e));
	p = {10, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,12) = " + std::to_string(e));
	p = {10, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,13) = " + std::to_string(e));
	p = {10, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,14) = " + std::to_string(e));
	p = {10, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,15) = " + std::to_string(e));
	p = {10, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,16) = " + std::to_string(e));
	p = {10, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(10,17) = " + std::to_string(e));
	
	p = {11, 11};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(11,11) = " + std::to_string(e));
	p = {11, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(11,12) = " + std::to_string(e));
	p = {11, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(11,13) = " + std::to_string(e));
	p = {11, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(11,14) = " + std::to_string(e));
	p = {11, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(11,15) = " + std::to_string(e));
	p = {11, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(11,16) = " + std::to_string(e));
	p = {11, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(11,17) = " + std::to_string(e));
	
	p = {12, 12};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(12,12) = " + std::to_string(e));
	p = {12, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(12,13) = " + std::to_string(e));
	p = {12, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(12,14) = " + std::to_string(e));
	p = {12, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(12,15) = " + std::to_string(e));
	p = {12, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(12,16) = " + std::to_string(e));
	p = {12, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(12,17) = " + std::to_string(e));
	
	p = {13, 13};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(13,13) = " + std::to_string(e));
	p = {13, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(13,14) = " + std::to_string(e));
	p = {13, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(13,15) = " + std::to_string(e));
	p = {13, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(13,16) = " + std::to_string(e));
	p = {13, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(13,17) = " + std::to_string(e));
	
	p = {14, 14};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(14,14) = " + std::to_string(e));
	p = {14, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(14,15) = " + std::to_string(e));
	p = {14, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(14,16) = " + std::to_string(e));
	p = {14, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(14,17) = " + std::to_string(e));
	
	p = {15, 15};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(15,15) = " + std::to_string(e));
	p = {15, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(15,16) = " + std::to_string(e));
	p = {15, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(15,17) = " + std::to_string(e));
	
	p = {16, 16};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(16,16) = " + std::to_string(e));
	p = {16, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(16,17) = " + std::to_string(e));
	
	p = {17, 17};
	e = zuker.wij.get(p);
	BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.00f), "W(17,17) = " + std::to_string(e));
	
}



void test_sij(Zuker &zuker, size_t n)
{
	Pair p = Pair();
	for(unsigned int x = 0;  x < (unsigned int) n; x++)
	{
		for(unsigned int y = x; y < (unsigned int) n; y++)
		{
			p = {x, y};
			BOOST_CHECK(zuker.sij.get(p) == nullptr);
		}
	}
}



/**
 * @date 2016-01-22
 */
void test_tij(Zuker &zuker)
{
	Pair p = Pair();
	traceback_jump2 t;
	
	p = {0, 0};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {0, 1};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {0, 2};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {0, 3};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {0, 4};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 4);
	p = {0, 5};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 5);
	p = {0, 6};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 6);
	p = {0, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 7);
	p = {0, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 8);
	p = {0, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {0, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {0, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {0, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 2);
	p = {0, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 12);
	p = {0, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {0, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {0, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {0, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 1);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	
	p = {1, 1};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {1, 2};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {1, 3};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {1, 4};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {1, 5};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 5);
	p = {1, 6};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 6);
	p = {1, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 7);
	p = {1, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 8);
	p = {1, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {1, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {1, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {1, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 3);
	p = {1, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 12);
	p = {1, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {1, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {1, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {1, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 2);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	
	p = {2, 2};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {2, 3};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {2, 4};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {2, 5};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {2, 6};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 6);
	p = {2, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 7);
	p = {2, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 8);
	p = {2, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {2, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {2, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {2, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 4);
	p = {2, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {2, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 4);
	p = {2, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 3);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {2, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 4);
	p = {2, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 4);
	
	p = {3, 3};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {3, 4};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {3, 5};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {3, 6};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {3, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 7);
	p = {3, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 8);
	p = {3, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {3, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {3, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {3, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 12);
	p = {3, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {3, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {3, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {3, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	p = {3, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 4);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {4, 4};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {4, 5};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {4, 6};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {4, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {4, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 8);
	p = {4, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {4, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {4, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {4, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 12);
	p = {4, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {4, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {4, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {4, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	p = {4, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 5);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {5, 5};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {5, 6};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {5, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {5, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {5, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {5, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {5, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {5, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 12);
	p = {5, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {5, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {5, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {5, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	p = {5, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 6);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {6, 6};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {6, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {6, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {6, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {6, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 7);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {6, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 7);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {6, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 8);
	BOOST_CHECK_EQUAL(t.target.second, 8);
	p = {6, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 7);
	BOOST_CHECK_EQUAL(t.target.second, 12);
	p = {6, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 7);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {6, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 7);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {6, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 7);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {6, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 7);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	
	p = {7, 7};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {7, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {7, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {7, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {7, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 8);
	BOOST_CHECK_EQUAL(t.target.second, 11);
	p = {7, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 9);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {7, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 9);
	BOOST_CHECK_EQUAL(t.target.second, 9);
	p = {7, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 8);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {7, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 8);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {7, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 8);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {7, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 8);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	
	p = {8, 8};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {8, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {8, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {8, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {8, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 9);
	BOOST_CHECK_EQUAL(t.target.second, 12);
	p = {8, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {8, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {8, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {8, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	p = {8, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, true);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 10);
	
	p = {9, 9};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {9, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {9, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {9, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {9, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 13);
	p = {9, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {9, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {9, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	p = {9, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 10);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {10, 10};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {10, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {10, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {10, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {10, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 11);
	BOOST_CHECK_EQUAL(t.target.second, 14);
	p = {10, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 11);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {10, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 11);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	p = {10, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 11);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {11, 11};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {11, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {11, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {11, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {11, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 12);
	BOOST_CHECK_EQUAL(t.target.second, 15);
	p = {11, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 12);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	p = {11, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 12);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {12, 12};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {12, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {12, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {12, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {12, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 13);
	BOOST_CHECK_EQUAL(t.target.second, 16);
	p = {12, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 13);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {13, 13};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {13, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {13, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {13, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {13, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, 14);
	BOOST_CHECK_EQUAL(t.target.second, 17);
	
	p = {14, 14};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {14, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {14, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {14, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	
	p = {15, 15};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {15, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {15, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	
	p = {16, 16};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	p = {16, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
	
	p = {17, 17};
	t = zuker.tij.get(p);
	BOOST_CHECK_EQUAL(t.store_pair, false);
	BOOST_CHECK_EQUAL(t.target.first, UNBOUND);
	//BOOST_CHECK_EQUAL(t.target.second, UNBOUND);
}



/**
 * @brief tests Bulge loop prediction
 *
 * @date 2016-01-22
 *
 * @test
 *
 * @section DESCRIPTION
 * Predicts folding and the ScoringMatrix contents for 1000 times running
 * the prediction.
 *
 * Initially, when the pathmatrix_corrected_from was a boolean container
 * some weird stuff hapenned.
 */
BOOST_AUTO_TEST_CASE(Test_bulge_loop)
{
	for(unsigned char i = 0 ; i < 10; i++)
	{
		// Initialize variables
		Sequence sequence = Sequence("GGGAAAGGGAAACCCCCC");
		std::string true_structure = "((....(((....)))))";
		
		Settings settings = Settings(0, nullptr, sequence);
		ReadData thermodynamics = ReadData();
		
		// Predict structure
		Zuker zuker = Zuker(settings, sequence, thermodynamics);
		zuker.energy();
		
		// The following code can be used to re-generate the funnctions above:
		/*
		//test_vij
		for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
		{
			for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
			{
				Pair p = Pair(x,y);
				printf("	p = {%i, %i};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, %.2ff),\"V(%i,%i) = \" + std::to_string(e) );\n", x,  y, zuker.vij.get(p), x, y);
			}
			printf("\n");
		}
		exit(1);
		
		//test_wij
		for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
		{
			for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
			{
				Pair p = Pair(x,y);
				printf("	p = {%i, %i};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, %.2ff),\"W(%i,%i) = \" + std::to_string(e) );\n", x,  y, zuker.wij.get(p), x, y);
			}
			printf("\n");
		}
		exit(1);
		//test_tij
		for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
		{
			for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
			{
				Pair p = Pair(x,y);
				traceback_jump2 t = zuker.tij.get(p);
				///todo - if t.target.first == UNBOUND  - do not check for t.target.second
				printf("	p = {%i, %i};t = zuker.tij.get(p);BOOST_CHECK_EQUAL(t.store_pair, %s);BOOST_CHECK_EQUAL(t.target.first, %i);BOOST_CHECK_EQUAL(t.target.second, %i);\n", x,  y, t.store_pair ? "true" : "false", t.target.first, t.target.second);
			}
			printf("\n");
		}
		exit(1);
		*/
		
		
		test_vij(zuker);
		test_wij(zuker);
		test_sij(zuker, sequence.size());
		test_tij(zuker);
		
		zuker.traceback();
		// Obtain and compare results
		std::string predicted_structure;
		zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure); ///@todo unsigned int -> size_t
		
		bool valid_structure = predicted_structure.compare(true_structure) == 0;
		BOOST_CHECK_MESSAGE(valid_structure, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	}
}

BOOST_AUTO_TEST_SUITE_END()
