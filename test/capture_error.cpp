/**
 * @file test/main_test.cpp
 *
 * @date 2015-07-23
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
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"
#include "ReadSegments.hpp"

#include "ScoringTree.hpp"

#include "ScoringMatrix.hpp"
#include "Settings.hpp"
#include "DotBracket.hpp"
#include "ReadData.hpp"

#include "Zuker.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)


bool AreSame(double a, double b)
{
    return fabs(a - b) < 0.0001;
}

void test_pij(Zuker &zuker)
{
	Pair p = Pair();
	p = {0,0};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {1,1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {2,2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {3,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {4,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {5,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {6,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {7,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {8,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {9,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {10,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {11,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {12,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {13,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {14,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {15,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {16,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {16,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {17,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
}

void test_qij(Zuker &zuker)
{
	Pair p = Pair();
	p = {0,0};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {1,1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {2,2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {3,3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {4,4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {5,5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {6,6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {7,7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {8,8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);

	p = {9,9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {10,10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {11,11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {12,12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {13,13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {14,14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {15,15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {16,16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {16,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);

	p = {17,17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);


}

void test_vij(Zuker &zuker)
{
	Pair p = Pair();
	
	p = {0,0};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {0,1};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {0,2};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {0,3};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {0,4};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,5};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,6};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {0,13};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {0,14};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {0,15};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {0,16};BOOST_CHECK(AreSame(zuker.vij.get(p), -9.900000));
	p = {0,17};BOOST_CHECK(AreSame(zuker.vij.get(p), -13.200000));

	p = {1,1};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {1,2};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {1,3};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {1,4};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {1,5};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,6};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {1,13};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {1,14};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {1,15};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {1,16};BOOST_CHECK(AreSame(zuker.vij.get(p), -9.900000));
	p = {1,17};BOOST_CHECK(AreSame(zuker.vij.get(p), -9.900000));

	p = {2,2};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {2,3};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {2,4};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {2,5};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {2,6};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {2,14};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {2,15};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {2,16};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {2,17};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));

	p = {3,3};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {3,4};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {3,5};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {3,6};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {3,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {3,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {4,4};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {4,5};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {4,6};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {4,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {4,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {4,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {5,5};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {5,6};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {5,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {5,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {5,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {5,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {6,6};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {6,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {6,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {6,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {6,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {6,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {6,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {6,13};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {6,14};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {6,15};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {6,16};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));
	p = {6,17};BOOST_CHECK(AreSame(zuker.vij.get(p), -6.600000));

	p = {7,7};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {7,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {7,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {7,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {7,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {7,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {7,13};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {7,14};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {7,15};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {7,16};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));
	p = {7,17};BOOST_CHECK(AreSame(zuker.vij.get(p), -3.300000));

	p = {8,8};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {8,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {8,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {8,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {8,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {8,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {8,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {8,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {8,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {8,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {9,9};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {9,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {9,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {9,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {9,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {9,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {9,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {9,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {9,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {10,10};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {10,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {10,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {10,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {10,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {10,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {10,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {10,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {11,11};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {11,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {11,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {11,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {11,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {11,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {11,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {12,12};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {12,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {12,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {12,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {12,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));
	p = {12,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {13,13};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {13,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {13,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {13,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {13,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 0.000000));

	p = {14,14};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {14,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {14,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {14,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));

	p = {15,15};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {15,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {15,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));

	p = {16,16};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));
	p = {16,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));

	p = {17,17};BOOST_CHECK(AreSame(zuker.vij.get(p), 999999.000000));

}

void test_wij(Zuker &zuker)
{
	Pair p = Pair();
	
	p = {0,0};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,1};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,2};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,3};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,4};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,5};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,6};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {0,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {0,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {0,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {0,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -9.900000));
	p = {0,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -13.200000));

	p = {1,1};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,2};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,3};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,4};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,5};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,6};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {1,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {1,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {1,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {1,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -9.900000));
	p = {1,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -9.900000));

	p = {2,2};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,3};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,4};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,5};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,6};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {2,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {2,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {2,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {2,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {2,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));

	p = {3,3};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,4};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,5};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,6};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {3,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {3,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {3,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {3,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {3,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));

	p = {4,4};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,5};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,6};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {4,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {4,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {4,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {4,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {4,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));

	p = {5,5};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,6};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {5,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {5,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {5,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {5,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {5,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));

	p = {6,6};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {6,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {6,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {6,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {6,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {6,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {6,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {6,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {6,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {6,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {6,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));
	p = {6,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -6.600000));

	p = {7,7};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {7,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {7,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {7,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {7,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {7,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {7,13};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {7,14};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {7,15};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {7,16};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));
	p = {7,17};BOOST_CHECK(AreSame(zuker.wij.get(p), -3.300000));

	p = {8,8};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,13};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,14};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {8,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {9,9};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,13};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,14};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {9,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {10,10};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {10,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {10,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {10,13};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {10,14};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {10,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {10,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {10,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {11,11};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {11,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {11,13};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {11,14};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {11,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {11,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {11,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {12,12};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {12,13};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {12,14};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {12,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {12,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {12,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {13,13};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {13,14};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {13,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {13,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {13,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {14,14};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {14,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {14,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {14,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {15,15};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {15,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {15,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {16,16};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
	p = {16,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));

	p = {17,17};BOOST_CHECK(AreSame(zuker.wij.get(p), 0.000000));
}

void test_pathmatrix_corrected_from(Zuker &zuker)
{
	Pair p = Pair();
	
	p = {0,0};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,1};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,2};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,3};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,4};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,5};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,6};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);//true);
	p = {0,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {0,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {0,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {0,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {0,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {0,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {0,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);

	p = {1,1};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,2};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,3};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,4};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,5};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,6};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {1,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {1,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {1,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {1,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {1,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);

	p = {2,2};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,3};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,4};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,5};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,6};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {2,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {2,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {2,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {2,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);

	p = {3,3};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,4};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,5};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,6};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);//true);
	p = {3,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {3,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {4,4};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,5};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,6};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {4,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {5,5};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,6};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {5,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {6,6};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {6,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {6,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {6,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {6,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {6,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {6,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {6,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {6,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {6,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {6,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {6,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);

	p = {7,7};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {7,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {7,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {7,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {7,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {7,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {7,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {7,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {7,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {7,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {7,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);

	p = {8,8};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {8,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {8,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {8,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {8,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {8,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {8,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {8,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {8,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);
	p = {8,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);

	p = {9,9};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {9,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), true);

	p = {10,10};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {10,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {10,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {10,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {10,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {10,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {10,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {10,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {11,11};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {11,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {11,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {11,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {11,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {11,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {11,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {12,12};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {12,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {12,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {12,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {12,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {12,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {13,13};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {13,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {13,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {13,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {13,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {14,14};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {14,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {14,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {14,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {15,15};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {15,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {15,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {16,16};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
	p = {16,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);

	p = {17,17};BOOST_CHECK_EQUAL(zuker.pathmatrix_corrected_from.get(p), false);
}

void test_loopmatrix(Zuker &zuker)
{
	Pair p = Pair();
	Pair p2 = Pair();

	p = {0,0};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,1};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,2};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 1);
	p = {0,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 12);
	p = {0,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {0,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {0,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {0,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {1,1};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,2};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 2);
	p = {1,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 12);
	p = {1,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {1,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {1,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {1,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {2,2};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 3);BOOST_CHECK_EQUAL(p2.second, 3);
	p = {2,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 3);BOOST_CHECK_EQUAL(p2.second, 3);
	p = {2,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 3);BOOST_CHECK_EQUAL(p2.second, 3);
	p = {2,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 3);BOOST_CHECK_EQUAL(p2.second, 3);
	p = {2,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 3);BOOST_CHECK_EQUAL(p2.second, 3);
	p = {2,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 3);BOOST_CHECK_EQUAL(p2.second, 3);

	p = {3,3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {4,4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {5,5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {6,6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 7);
	p = {6,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 12);
	p = {6,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {6,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {6,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {6,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {7,7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 8);
	p = {7,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 12);
	p = {7,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {7,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {7,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {7,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {8,8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);
	p = {8,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);
	p = {8,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);
	p = {8,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);
	p = {8,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);
	p = {8,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);

	p = {9,9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {10,10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {11,11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {12,12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {13,13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {14,14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {14,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {14,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {14,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {15,15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {15,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {15,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {16,16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {16,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {17,17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
}

void test_nij2(Zuker &zuker)
{
	Pair p = Pair();

	p = {0,0};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,1};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,2};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,3};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,4};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,5};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,6};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {0,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {1,1};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,2};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,3};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,4};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,5};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,6};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {1,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {2,2};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,3};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,4};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,5};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,6};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {2,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {3,3};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,4};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,5};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,6};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {3,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {4,4};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,5};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,6};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {4,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {5,5};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,6};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {5,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {6,6};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {6,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {7,7};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {7,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {8,8};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {8,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {9,9};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {9,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {10,10};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {10,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {10,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {10,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {10,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {10,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {10,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {10,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {11,11};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {11,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {11,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {11,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {11,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {11,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {11,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {12,12};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {12,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {12,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {12,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {12,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {12,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {13,13};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {13,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {13,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {13,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {13,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {14,14};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {14,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {14,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {14,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {15,15};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {15,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {15,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {16,16};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
	p = {16,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);

	p = {17,17};BOOST_CHECK(zuker.nij2.get(p) == nullptr);
}

/**
 * @brief tests Bulge loop prediction
 *
 * @test
 *
 * @section DESCRIPTION
 * Predicts folding of the following bulge-loop structure:
 * <PRE>
 * 5') GGGAAAGGG A
 *     \\\   ///  A
 * 3') CCC CCC A
 * </PRE>
 */
BOOST_AUTO_TEST_CASE(Test_bulge_loop)
{
	int f =  0;
	for(int i = 0 ; i < 1000; i++)
	{
		// Initialize variables
		Sequence sequence = Sequence("GGGAAAGGGAAACCCCCC");
		std::string true_structure = "(((...(((...))))))";
		
		Settings settings = Settings(0, nullptr, sequence);
		ReadData thermodynamics = ReadData();
		
		// Predict structure
		Zuker zuker = Zuker(settings, sequence, thermodynamics);
		zuker.energy();
		zuker.traceback();
		
		// Obtain and compare results
		std::string predicted_structure;
		zuker.dot_bracket.format(sequence.size() , predicted_structure);
		
		bool valid_structure = predicted_structure.compare(true_structure) == 0;
		BOOST_CHECK_MESSAGE(valid_structure, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
		
		if(valid_structure)
		{
			f++;
			test_pij(zuker);
			test_qij(zuker);
			test_vij(zuker);
			test_wij(zuker);
			test_nij2(zuker);
			test_pathmatrix_corrected_from(zuker);
			test_loopmatrix(zuker);
			
			
			//pathmatrix_corrected_from
			/*
			for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
			{
				for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
				{
					Pair p = Pair(x,y);
					Segment *n = zuker.nij2.get(p);
					if(n == nullptr || n == NULL)
					{
						printf("	p = {%i,%i};BOOST_CHECK_EQUAL(zuker.nij2.get(p), nullptr);\n", x,  y);
					}
					else
					{
						printf("	p = {%i,%i};BOOST_CHECK(zuker.nij2.get(p) != nullptr);\n", x,  y);
					}
				}
				printf("\n");
			}
			*/
		}
	}
}



BOOST_AUTO_TEST_SUITE_END()
