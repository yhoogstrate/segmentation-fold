/**
 * @file test/Zuker_traceback_test.cpp
 *
 * @date 2015-12-01
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


bool IdenticalFloats(double a, double b)
{
	return fabs(a - b) < 0.0001;
}

void test_pij(Zuker &zuker)
{
	Pair p = Pair();
	
	p = {0, 0};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {1, 1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {2, 2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {3, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {4, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {5, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {6, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {7, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {8, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {9, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {10, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {11, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {12, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {13, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {14, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {15, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {16, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {16, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {17, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
}

void test_qij(Zuker &zuker)
{
	Pair p = Pair();
	p = {0, 0};
	BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {0, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {0, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {1, 1};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {1, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {1, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {2, 2};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {2, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {2, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {3, 3};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {3, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {4, 4};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {4, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {5, 5};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {5, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {6, 6};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {6, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {6, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {7, 7};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {7, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {7, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {8, 8};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {8, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	p = {8, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -1);
	
	p = {9, 9};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {9, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {10, 10};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {10, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {11, 11};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {11, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {12, 12};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {12, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {13, 13};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {13, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {14, 14};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {14, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {15, 15};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {15, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {16, 16};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	p = {16, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
	
	p = {17, 17};BOOST_CHECK_EQUAL(zuker.qij.get(p), -2);
}

void test_vij(Zuker &zuker)
{
	Pair p = Pair();
	float e;
	
	p = {0, 0};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(0,0) = " + std::to_string(e) );
	p = {0, 1};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(0,1) = " + std::to_string(e) );
	p = {0, 2};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(0,2) = " + std::to_string(e) );
	p = {0, 3};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(0,3) = " + std::to_string(e) );
	p = {0, 4};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,4) = " + std::to_string(e) );
	p = {0, 5};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,5) = " + std::to_string(e) );
	p = {0, 6};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,6) = " + std::to_string(e) );
	p = {0, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,7) = " + std::to_string(e) );
	p = {0, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,8) = " + std::to_string(e) );
	p = {0, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,9) = " + std::to_string(e) );
	p = {0, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,10) = " + std::to_string(e) );
	p = {0, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,11) = " + std::to_string(e) );
	p = {0, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(0,12) = " + std::to_string(e) );
	p = {0, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(0,13) = " + std::to_string(e) );
	p = {0, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(0,14) = " + std::to_string(e) );
	p = {0, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(0,15) = " + std::to_string(e) );
	p = {0, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.900000),"V(0,16) = " + std::to_string(e) );
	p = {0, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.900000),"V(0,17) = " + std::to_string(e) );

	p = {1, 1};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(1,1) = " + std::to_string(e) );
	p = {1, 2};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(1,2) = " + std::to_string(e) );
	p = {1, 3};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(1,3) = " + std::to_string(e) );
	p = {1, 4};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(1,4) = " + std::to_string(e) );
	p = {1, 5};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,5) = " + std::to_string(e) );
	p = {1, 6};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,6) = " + std::to_string(e) );
	p = {1, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,7) = " + std::to_string(e) );
	p = {1, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,8) = " + std::to_string(e) );
	p = {1, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,9) = " + std::to_string(e) );
	p = {1, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,10) = " + std::to_string(e) );
	p = {1, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,11) = " + std::to_string(e) );
	p = {1, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(1,12) = " + std::to_string(e) );
	p = {1, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(1,13) = " + std::to_string(e) );
	p = {1, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(1,14) = " + std::to_string(e) );
	p = {1, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(1,15) = " + std::to_string(e) );
	p = {1, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(1,16) = " + std::to_string(e) );
	p = {1, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.900000),"V(1,17) = " + std::to_string(e) );

	p = {2, 2};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(2,2) = " + std::to_string(e) );
	p = {2, 3};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(2,3) = " + std::to_string(e) );
	p = {2, 4};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(2,4) = " + std::to_string(e) );
	p = {2, 5};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(2,5) = " + std::to_string(e) );
	p = {2, 6};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,6) = " + std::to_string(e) );
	p = {2, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,7) = " + std::to_string(e) );
	p = {2, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,8) = " + std::to_string(e) );
	p = {2, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,9) = " + std::to_string(e) );
	p = {2, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,10) = " + std::to_string(e) );
	p = {2, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,11) = " + std::to_string(e) );
	p = {2, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,12) = " + std::to_string(e) );
	p = {2, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(2,13) = " + std::to_string(e) );
	p = {2, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(2,14) = " + std::to_string(e) );
	p = {2, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(2,15) = " + std::to_string(e) );
	p = {2, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(2,16) = " + std::to_string(e) );
	p = {2, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(2,17) = " + std::to_string(e) );

	p = {3, 3};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(3,3) = " + std::to_string(e) );
	p = {3, 4};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(3,4) = " + std::to_string(e) );
	p = {3, 5};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(3,5) = " + std::to_string(e) );
	p = {3, 6};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(3,6) = " + std::to_string(e) );
	p = {3, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,7) = " + std::to_string(e) );
	p = {3, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,8) = " + std::to_string(e) );
	p = {3, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,9) = " + std::to_string(e) );
	p = {3, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,10) = " + std::to_string(e) );
	p = {3, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,11) = " + std::to_string(e) );
	p = {3, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,12) = " + std::to_string(e) );
	p = {3, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,13) = " + std::to_string(e) );
	p = {3, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,14) = " + std::to_string(e) );
	p = {3, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,15) = " + std::to_string(e) );
	p = {3, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,16) = " + std::to_string(e) );
	p = {3, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(3,17) = " + std::to_string(e) );

	p = {4, 4};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(4,4) = " + std::to_string(e) );
	p = {4, 5};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(4,5) = " + std::to_string(e) );
	p = {4, 6};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(4,6) = " + std::to_string(e) );
	p = {4, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(4,7) = " + std::to_string(e) );
	p = {4, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,8) = " + std::to_string(e) );
	p = {4, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,9) = " + std::to_string(e) );
	p = {4, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,10) = " + std::to_string(e) );
	p = {4, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,11) = " + std::to_string(e) );
	p = {4, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,12) = " + std::to_string(e) );
	p = {4, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,13) = " + std::to_string(e) );
	p = {4, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,14) = " + std::to_string(e) );
	p = {4, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,15) = " + std::to_string(e) );
	p = {4, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,16) = " + std::to_string(e) );
	p = {4, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(4,17) = " + std::to_string(e) );

	p = {5, 5};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(5,5) = " + std::to_string(e) );
	p = {5, 6};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(5,6) = " + std::to_string(e) );
	p = {5, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(5,7) = " + std::to_string(e) );
	p = {5, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(5,8) = " + std::to_string(e) );
	p = {5, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,9) = " + std::to_string(e) );
	p = {5, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,10) = " + std::to_string(e) );
	p = {5, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,11) = " + std::to_string(e) );
	p = {5, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,12) = " + std::to_string(e) );
	p = {5, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,13) = " + std::to_string(e) );
	p = {5, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,14) = " + std::to_string(e) );
	p = {5, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,15) = " + std::to_string(e) );
	p = {5, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,16) = " + std::to_string(e) );
	p = {5, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(5,17) = " + std::to_string(e) );

	p = {6, 6};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(6,6) = " + std::to_string(e) );
	p = {6, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(6,7) = " + std::to_string(e) );
	p = {6, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(6,8) = " + std::to_string(e) );
	p = {6, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(6,9) = " + std::to_string(e) );
	p = {6, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(6,10) = " + std::to_string(e) );
	p = {6, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(6,11) = " + std::to_string(e) );
	p = {6, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(6,12) = " + std::to_string(e) );
	p = {6, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(6,13) = " + std::to_string(e) );
	p = {6, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(6,14) = " + std::to_string(e) );
	p = {6, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(6,15) = " + std::to_string(e) );
	p = {6, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(6,16) = " + std::to_string(e) );
	p = {6, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"V(6,17) = " + std::to_string(e) );

	p = {7, 7};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(7,7) = " + std::to_string(e) );
	p = {7, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(7,8) = " + std::to_string(e) );
	p = {7, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(7,9) = " + std::to_string(e) );
	p = {7, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(7,10) = " + std::to_string(e) );
	p = {7, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(7,11) = " + std::to_string(e) );
	p = {7, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(7,12) = " + std::to_string(e) );
	p = {7, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(7,13) = " + std::to_string(e) );
	p = {7, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(7,14) = " + std::to_string(e) );
	p = {7, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(7,15) = " + std::to_string(e) );
	p = {7, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(7,16) = " + std::to_string(e) );
	p = {7, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"V(7,17) = " + std::to_string(e) );

	p = {8, 8};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(8,8) = " + std::to_string(e) );
	p = {8, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(8,9) = " + std::to_string(e) );
	p = {8, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(8,10) = " + std::to_string(e) );
	p = {8, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(8,11) = " + std::to_string(e) );
	p = {8, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 4.600000),"V(8,12) = " + std::to_string(e) );
	p = {8, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(8,13) = " + std::to_string(e) );
	p = {8, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(8,14) = " + std::to_string(e) );
	p = {8, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(8,15) = " + std::to_string(e) );
	p = {8, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(8,16) = " + std::to_string(e) );
	p = {8, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(8,17) = " + std::to_string(e) );

	p = {9, 9};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(9,9) = " + std::to_string(e) );
	p = {9, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(9,10) = " + std::to_string(e) );
	p = {9, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(9,11) = " + std::to_string(e) );
	p = {9, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(9,12) = " + std::to_string(e) );
	p = {9, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(9,13) = " + std::to_string(e) );
	p = {9, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(9,14) = " + std::to_string(e) );
	p = {9, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(9,15) = " + std::to_string(e) );
	p = {9, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(9,16) = " + std::to_string(e) );
	p = {9, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(9,17) = " + std::to_string(e) );

	p = {10, 10};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(10,10) = " + std::to_string(e) );
	p = {10, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(10,11) = " + std::to_string(e) );
	p = {10, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(10,12) = " + std::to_string(e) );
	p = {10, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(10,13) = " + std::to_string(e) );
	p = {10, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(10,14) = " + std::to_string(e) );
	p = {10, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(10,15) = " + std::to_string(e) );
	p = {10, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(10,16) = " + std::to_string(e) );
	p = {10, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(10,17) = " + std::to_string(e) );

	p = {11, 11};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(11,11) = " + std::to_string(e) );
	p = {11, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(11,12) = " + std::to_string(e) );
	p = {11, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(11,13) = " + std::to_string(e) );
	p = {11, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(11,14) = " + std::to_string(e) );
	p = {11, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(11,15) = " + std::to_string(e) );
	p = {11, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(11,16) = " + std::to_string(e) );
	p = {11, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(11,17) = " + std::to_string(e) );

	p = {12, 12};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(12,12) = " + std::to_string(e) );
	p = {12, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(12,13) = " + std::to_string(e) );
	p = {12, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(12,14) = " + std::to_string(e) );
	p = {12, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(12,15) = " + std::to_string(e) );
	p = {12, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(12,16) = " + std::to_string(e) );
	p = {12, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(12,17) = " + std::to_string(e) );

	p = {13, 13};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(13,13) = " + std::to_string(e) );
	p = {13, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(13,14) = " + std::to_string(e) );
	p = {13, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(13,15) = " + std::to_string(e) );
	p = {13, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(13,16) = " + std::to_string(e) );
	p = {13, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"V(13,17) = " + std::to_string(e) );

	p = {14, 14};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(14,14) = " + std::to_string(e) );
	p = {14, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(14,15) = " + std::to_string(e) );
	p = {14, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(14,16) = " + std::to_string(e) );
	p = {14, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(14,17) = " + std::to_string(e) );

	p = {15, 15};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(15,15) = " + std::to_string(e) );
	p = {15, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(15,16) = " + std::to_string(e) );
	p = {15, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(15,17) = " + std::to_string(e) );

	p = {16, 16};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(16,16) = " + std::to_string(e) );
	p = {16, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(16,17) = " + std::to_string(e) );

	p = {17, 17};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 999999.000000),"V(17,17) = " + std::to_string(e) );
}

void test_wij(Zuker &zuker)
{
	Pair p = Pair();
	float e;
	
	p = {0, 0};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,0) = " + std::to_string(e) );
	p = {0, 1};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,1) = " + std::to_string(e) );
	p = {0, 2};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,2) = " + std::to_string(e) );
	p = {0, 3};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,3) = " + std::to_string(e) );
	p = {0, 4};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,4) = " + std::to_string(e) );
	p = {0, 5};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,5) = " + std::to_string(e) );
	p = {0, 6};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,6) = " + std::to_string(e) );
	p = {0, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,7) = " + std::to_string(e) );
	p = {0, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,8) = " + std::to_string(e) );
	p = {0, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,9) = " + std::to_string(e) );
	p = {0, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,10) = " + std::to_string(e) );
	p = {0, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,11) = " + std::to_string(e) );
	p = {0, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(0,12) = " + std::to_string(e) );
	p = {0, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(0,13) = " + std::to_string(e) );
	p = {0, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(0,14) = " + std::to_string(e) );
	p = {0, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(0,15) = " + std::to_string(e) );
	p = {0, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.900000),"W(0,16) = " + std::to_string(e) );
	p = {0, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.900000),"W(0,17) = " + std::to_string(e) );

	p = {1, 1};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,1) = " + std::to_string(e) );
	p = {1, 2};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,2) = " + std::to_string(e) );
	p = {1, 3};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,3) = " + std::to_string(e) );
	p = {1, 4};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,4) = " + std::to_string(e) );
	p = {1, 5};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,5) = " + std::to_string(e) );
	p = {1, 6};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,6) = " + std::to_string(e) );
	p = {1, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,7) = " + std::to_string(e) );
	p = {1, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,8) = " + std::to_string(e) );
	p = {1, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,9) = " + std::to_string(e) );
	p = {1, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,10) = " + std::to_string(e) );
	p = {1, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,11) = " + std::to_string(e) );
	p = {1, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(1,12) = " + std::to_string(e) );
	p = {1, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(1,13) = " + std::to_string(e) );
	p = {1, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(1,14) = " + std::to_string(e) );
	p = {1, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(1,15) = " + std::to_string(e) );
	p = {1, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(1,16) = " + std::to_string(e) );
	p = {1, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -9.900000),"W(1,17) = " + std::to_string(e) );

	p = {2, 2};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,2) = " + std::to_string(e) );
	p = {2, 3};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,3) = " + std::to_string(e) );
	p = {2, 4};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,4) = " + std::to_string(e) );
	p = {2, 5};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,5) = " + std::to_string(e) );
	p = {2, 6};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,6) = " + std::to_string(e) );
	p = {2, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,7) = " + std::to_string(e) );
	p = {2, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,8) = " + std::to_string(e) );
	p = {2, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,9) = " + std::to_string(e) );
	p = {2, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,10) = " + std::to_string(e) );
	p = {2, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,11) = " + std::to_string(e) );
	p = {2, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(2,12) = " + std::to_string(e) );
	p = {2, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(2,13) = " + std::to_string(e) );
	p = {2, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(2,14) = " + std::to_string(e) );
	p = {2, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(2,15) = " + std::to_string(e) );
	p = {2, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(2,16) = " + std::to_string(e) );
	p = {2, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(2,17) = " + std::to_string(e) );

	p = {3, 3};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,3) = " + std::to_string(e) );
	p = {3, 4};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,4) = " + std::to_string(e) );
	p = {3, 5};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,5) = " + std::to_string(e) );
	p = {3, 6};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,6) = " + std::to_string(e) );
	p = {3, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,7) = " + std::to_string(e) );
	p = {3, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,8) = " + std::to_string(e) );
	p = {3, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,9) = " + std::to_string(e) );
	p = {3, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,10) = " + std::to_string(e) );
	p = {3, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,11) = " + std::to_string(e) );
	p = {3, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(3,12) = " + std::to_string(e) );
	p = {3, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(3,13) = " + std::to_string(e) );
	p = {3, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(3,14) = " + std::to_string(e) );
	p = {3, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(3,15) = " + std::to_string(e) );
	p = {3, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(3,16) = " + std::to_string(e) );
	p = {3, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(3,17) = " + std::to_string(e) );

	p = {4, 4};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,4) = " + std::to_string(e) );
	p = {4, 5};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,5) = " + std::to_string(e) );
	p = {4, 6};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,6) = " + std::to_string(e) );
	p = {4, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,7) = " + std::to_string(e) );
	p = {4, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,8) = " + std::to_string(e) );
	p = {4, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,9) = " + std::to_string(e) );
	p = {4, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,10) = " + std::to_string(e) );
	p = {4, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,11) = " + std::to_string(e) );
	p = {4, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(4,12) = " + std::to_string(e) );
	p = {4, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(4,13) = " + std::to_string(e) );
	p = {4, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(4,14) = " + std::to_string(e) );
	p = {4, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(4,15) = " + std::to_string(e) );
	p = {4, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(4,16) = " + std::to_string(e) );
	p = {4, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(4,17) = " + std::to_string(e) );

	p = {5, 5};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,5) = " + std::to_string(e) );
	p = {5, 6};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,6) = " + std::to_string(e) );
	p = {5, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,7) = " + std::to_string(e) );
	p = {5, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,8) = " + std::to_string(e) );
	p = {5, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,9) = " + std::to_string(e) );
	p = {5, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,10) = " + std::to_string(e) );
	p = {5, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,11) = " + std::to_string(e) );
	p = {5, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,12) = " + std::to_string(e) );
	p = {5, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(5,13) = " + std::to_string(e) );
	p = {5, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(5,14) = " + std::to_string(e) );
	p = {5, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(5,15) = " + std::to_string(e) );
	p = {5, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(5,16) = " + std::to_string(e) );
	p = {5, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(5,17) = " + std::to_string(e) );

	p = {6, 6};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(6,6) = " + std::to_string(e) );
	p = {6, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(6,7) = " + std::to_string(e) );
	p = {6, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(6,8) = " + std::to_string(e) );
	p = {6, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(6,9) = " + std::to_string(e) );
	p = {6, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(6,10) = " + std::to_string(e) );
	p = {6, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(6,11) = " + std::to_string(e) );
	p = {6, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(6,12) = " + std::to_string(e) );
	p = {6, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(6,13) = " + std::to_string(e) );
	p = {6, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(6,14) = " + std::to_string(e) );
	p = {6, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(6,15) = " + std::to_string(e) );
	p = {6, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(6,16) = " + std::to_string(e) );
	p = {6, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -6.600000),"W(6,17) = " + std::to_string(e) );

	p = {7, 7};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(7,7) = " + std::to_string(e) );
	p = {7, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(7,8) = " + std::to_string(e) );
	p = {7, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(7,9) = " + std::to_string(e) );
	p = {7, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(7,10) = " + std::to_string(e) );
	p = {7, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(7,11) = " + std::to_string(e) );
	p = {7, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(7,12) = " + std::to_string(e) );
	p = {7, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(7,13) = " + std::to_string(e) );
	p = {7, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(7,14) = " + std::to_string(e) );
	p = {7, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(7,15) = " + std::to_string(e) );
	p = {7, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(7,16) = " + std::to_string(e) );
	p = {7, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, -3.300000),"W(7,17) = " + std::to_string(e) );

	p = {8, 8};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,8) = " + std::to_string(e) );
	p = {8, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,9) = " + std::to_string(e) );
	p = {8, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,10) = " + std::to_string(e) );
	p = {8, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,11) = " + std::to_string(e) );
	p = {8, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,12) = " + std::to_string(e) );
	p = {8, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,13) = " + std::to_string(e) );
	p = {8, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,14) = " + std::to_string(e) );
	p = {8, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,15) = " + std::to_string(e) );
	p = {8, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,16) = " + std::to_string(e) );
	p = {8, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(8,17) = " + std::to_string(e) );

	p = {9, 9};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,9) = " + std::to_string(e) );
	p = {9, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,10) = " + std::to_string(e) );
	p = {9, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,11) = " + std::to_string(e) );
	p = {9, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,12) = " + std::to_string(e) );
	p = {9, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,13) = " + std::to_string(e) );
	p = {9, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,14) = " + std::to_string(e) );
	p = {9, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,15) = " + std::to_string(e) );
	p = {9, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,16) = " + std::to_string(e) );
	p = {9, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(9,17) = " + std::to_string(e) );

	p = {10, 10};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,10) = " + std::to_string(e) );
	p = {10, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,11) = " + std::to_string(e) );
	p = {10, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,12) = " + std::to_string(e) );
	p = {10, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,13) = " + std::to_string(e) );
	p = {10, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,14) = " + std::to_string(e) );
	p = {10, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,15) = " + std::to_string(e) );
	p = {10, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,16) = " + std::to_string(e) );
	p = {10, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(10,17) = " + std::to_string(e) );

	p = {11, 11};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(11,11) = " + std::to_string(e) );
	p = {11, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(11,12) = " + std::to_string(e) );
	p = {11, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(11,13) = " + std::to_string(e) );
	p = {11, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(11,14) = " + std::to_string(e) );
	p = {11, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(11,15) = " + std::to_string(e) );
	p = {11, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(11,16) = " + std::to_string(e) );
	p = {11, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(11,17) = " + std::to_string(e) );

	p = {12, 12};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(12,12) = " + std::to_string(e) );
	p = {12, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(12,13) = " + std::to_string(e) );
	p = {12, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(12,14) = " + std::to_string(e) );
	p = {12, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(12,15) = " + std::to_string(e) );
	p = {12, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(12,16) = " + std::to_string(e) );
	p = {12, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(12,17) = " + std::to_string(e) );

	p = {13, 13};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(13,13) = " + std::to_string(e) );
	p = {13, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(13,14) = " + std::to_string(e) );
	p = {13, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(13,15) = " + std::to_string(e) );
	p = {13, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(13,16) = " + std::to_string(e) );
	p = {13, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(13,17) = " + std::to_string(e) );

	p = {14, 14};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(14,14) = " + std::to_string(e) );
	p = {14, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(14,15) = " + std::to_string(e) );
	p = {14, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(14,16) = " + std::to_string(e) );
	p = {14, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(14,17) = " + std::to_string(e) );

	p = {15, 15};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(15,15) = " + std::to_string(e) );
	p = {15, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(15,16) = " + std::to_string(e) );
	p = {15, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(15,17) = " + std::to_string(e) );

	p = {16, 16};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(16,16) = " + std::to_string(e) );
	p = {16, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(16,17) = " + std::to_string(e) );

	p = {17, 17};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, 0.000000),"W(17,17) = " + std::to_string(e) );
}

void test_pathmatrix_corrected_from(Zuker &zuker)
{
	Pair p = Pair();
	
	p = {0, 0};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,0) = true" );
	p = {0, 1};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,1) = true" );
	p = {0, 2};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,2) = true" );
	p = {0, 3};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,3) = true" );
	p = {0, 4};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,4) = true" );
	p = {0, 5};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,5) = true" );
	p = {0, 6};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,6) = true" );
	p = {0, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,7) = true" );
	p = {0, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,8) = true" );
	p = {0, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,9) = true" );
	p = {0, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,10) = true" );
	p = {0, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(0,11) = true" );
	p = {0, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(0,12) = false" );
	p = {0, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(0,13) = false" );
	p = {0, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(0,14) = false" );
	p = {0, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(0,15) = false" );
	p = {0, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(0,16) = false" );
	p = {0, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(0,17) = false" );

	p = {1, 1};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,1) = true" );
	p = {1, 2};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,2) = true" );
	p = {1, 3};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,3) = true" );
	p = {1, 4};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,4) = true" );
	p = {1, 5};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,5) = true" );
	p = {1, 6};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,6) = true" );
	p = {1, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,7) = true" );
	p = {1, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,8) = true" );
	p = {1, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,9) = true" );
	p = {1, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,10) = true" );
	p = {1, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(1,11) = true" );
	p = {1, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(1,12) = false" );
	p = {1, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(1,13) = false" );
	p = {1, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(1,14) = false" );
	p = {1, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(1,15) = false" );
	p = {1, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(1,16) = false" );
	p = {1, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(1,17) = false" );

	p = {2, 2};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,2) = true" );
	p = {2, 3};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,3) = true" );
	p = {2, 4};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,4) = true" );
	p = {2, 5};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,5) = true" );
	p = {2, 6};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,6) = true" );
	p = {2, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,7) = true" );
	p = {2, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,8) = true" );
	p = {2, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,9) = true" );
	p = {2, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,10) = true" );
	p = {2, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,11) = true" );
	p = {2, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(2,12) = false" );
	p = {2, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,13) = true" );
	p = {2, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(2,14) = false" );
	p = {2, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(2,15) = true" );
	p = {2, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(2,16) = false" );
	p = {2, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(2,17) = false" );

	p = {3, 3};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,3) = true" );
	p = {3, 4};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,4) = true" );
	p = {3, 5};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,5) = true" );
	p = {3, 6};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,6) = true" );
	p = {3, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,7) = true" );
	p = {3, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,8) = true" );
	p = {3, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,9) = true" );
	p = {3, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,10) = true" );
	p = {3, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,11) = true" );
	p = {3, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,12) = true" );
	p = {3, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,13) = true" );
	p = {3, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,14) = true" );
	p = {3, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,15) = true" );
	p = {3, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,16) = true" );
	p = {3, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(3,17) = true" );

	p = {4, 4};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,4) = true" );
	p = {4, 5};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,5) = true" );
	p = {4, 6};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,6) = true" );
	p = {4, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,7) = true" );
	p = {4, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,8) = true" );
	p = {4, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,9) = true" );
	p = {4, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,10) = true" );
	p = {4, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,11) = true" );
	p = {4, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,12) = true" );
	p = {4, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,13) = true" );
	p = {4, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,14) = true" );
	p = {4, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,15) = true" );
	p = {4, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,16) = true" );
	p = {4, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(4,17) = true" );

	p = {5, 5};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,5) = true" );
	p = {5, 6};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,6) = true" );
	p = {5, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,7) = true" );
	p = {5, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,8) = true" );
	p = {5, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,9) = true" );
	p = {5, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,10) = true" );
	p = {5, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,11) = true" );
	p = {5, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,12) = true" );
	p = {5, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,13) = true" );
	p = {5, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,14) = true" );
	p = {5, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,15) = true" );
	p = {5, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,16) = true" );
	p = {5, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(5,17) = true" );

	p = {6, 6};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(6,6) = true" );
	p = {6, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(6,7) = true" );
	p = {6, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(6,8) = true" );
	p = {6, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(6,9) = true" );
	p = {6, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(6,10) = true" );
	p = {6, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(6,11) = true" );
	p = {6, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(6,12) = false" );
	p = {6, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(6,13) = false" );
	p = {6, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(6,14) = false" );
	p = {6, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(6,15) = false" );
	p = {6, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(6,16) = false" );
	p = {6, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(6,17) = false" );

	p = {7, 7};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(7,7) = true" );
	p = {7, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(7,8) = true" );
	p = {7, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(7,9) = true" );
	p = {7, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(7,10) = true" );
	p = {7, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(7,11) = true" );
	p = {7, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(7,12) = false" );
	p = {7, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(7,13) = false" );
	p = {7, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(7,14) = false" );
	p = {7, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(7,15) = false" );
	p = {7, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(7,16) = false" );
	p = {7, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(7,17) = false" );

	p = {8, 8};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(8,8) = true" );
	p = {8, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(8,9) = true" );
	p = {8, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(8,10) = true" );
	p = {8, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(8,11) = true" );
	p = {8, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(8,12) = true" );
	p = {8, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(8,13) = false" );
	p = {8, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(8,14) = false" );
	p = {8, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(8,15) = false" );
	p = {8, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(8,16) = false" );
	p = {8, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,"path_matrix_corrected_from(8,17) = false" );

	p = {9, 9};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,9) = true" );
	p = {9, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,10) = true" );
	p = {9, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,11) = true" );
	p = {9, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,12) = true" );
	p = {9, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,13) = true" );
	p = {9, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,14) = true" );
	p = {9, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,15) = true" );
	p = {9, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,16) = true" );
	p = {9, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(9,17) = true" );

	p = {10, 10};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,10) = true" );
	p = {10, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,11) = true" );
	p = {10, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,12) = true" );
	p = {10, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,13) = true" );
	p = {10, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,14) = true" );
	p = {10, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,15) = true" );
	p = {10, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,16) = true" );
	p = {10, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(10,17) = true" );

	p = {11, 11};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(11,11) = true" );
	p = {11, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(11,12) = true" );
	p = {11, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(11,13) = true" );
	p = {11, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(11,14) = true" );
	p = {11, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(11,15) = true" );
	p = {11, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(11,16) = true" );
	p = {11, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(11,17) = true" );

	p = {12, 12};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(12,12) = true" );
	p = {12, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(12,13) = true" );
	p = {12, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(12,14) = true" );
	p = {12, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(12,15) = true" );
	p = {12, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(12,16) = true" );
	p = {12, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(12,17) = true" );

	p = {13, 13};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(13,13) = true" );
	p = {13, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(13,14) = true" );
	p = {13, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(13,15) = true" );
	p = {13, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(13,16) = true" );
	p = {13, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(13,17) = true" );

	p = {14, 14};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(14,14) = true" );
	p = {14, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(14,15) = true" );
	p = {14, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(14,16) = true" );
	p = {14, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(14,17) = true" );

	p = {15, 15};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(15,15) = true" );
	p = {15, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(15,16) = true" );
	p = {15, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(15,17) = true" );

	p = {16, 16};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(16,16) = true" );
	p = {16, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(16,17) = true" );

	p = {17, 17};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,"path_matrix_corrected_from(17,17) = true" );
}

void test_loopmatrix(Zuker &zuker)
{
	Pair p = Pair();
	Pair p2 = Pair();
	
	p = {0, 0};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 1};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 2};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {0, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 2);
	p = {0, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 12);
	p = {0, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {0, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {0, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {0, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 1);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {1, 1};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 2};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {1, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 3);BOOST_CHECK_EQUAL(p2.second, 3);
	p = {1, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 12);
	p = {1, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {1, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {1, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {1, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 2);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {2, 2};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {2, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 4);BOOST_CHECK_EQUAL(p2.second, 4);
	p = {2, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 4);BOOST_CHECK_EQUAL(p2.second, 4);
	p = {2, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 5);BOOST_CHECK_EQUAL(p2.second, 5);
	p = {2, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 4);BOOST_CHECK_EQUAL(p2.second, 4);
	p = {2, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 5);BOOST_CHECK_EQUAL(p2.second, 5);
	p = {2, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 5);BOOST_CHECK_EQUAL(p2.second, 5);

	p = {3, 3};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {3, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {4, 4};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {4, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {5, 5};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {5, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {6, 6};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {6, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 8);
	p = {6, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 12);
	p = {6, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {6, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {6, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {6, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 7);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {7, 7};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {7, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);
	p = {7, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 9);
	p = {7, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 13);
	p = {7, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 14);
	p = {7, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 15);
	p = {7, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 8);BOOST_CHECK_EQUAL(p2.second, 16);

	p = {8, 8};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {8, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 9);BOOST_CHECK_EQUAL(p2.second, 11);
	p = {8, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 10);BOOST_CHECK_EQUAL(p2.second, 10);
	p = {8, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 10);BOOST_CHECK_EQUAL(p2.second, 10);
	p = {8, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 10);BOOST_CHECK_EQUAL(p2.second, 10);
	p = {8, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 10);BOOST_CHECK_EQUAL(p2.second, 10);
	p = {8, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 10);BOOST_CHECK_EQUAL(p2.second, 10);

	p = {9, 9};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {9, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {10, 10};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {10, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {11, 11};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {11, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {12, 12};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {12, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {13, 13};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {13, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {14, 14};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {14, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {14, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {14, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {15, 15};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {15, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {15, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {16, 16};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
	p = {16, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);

	p = {17, 17};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, 0);BOOST_CHECK_EQUAL(p2.second, 0);
}

void test_nij2(Zuker &zuker, size_t n)
{
	Pair p = Pair();
	for(unsigned int x = 0;  x < (unsigned int) n; x++)
	{
		for(unsigned int y = x; y < (unsigned int) n; y++)
		{
			p = {x, y};
			BOOST_CHECK(zuker.nij2.get(p) == nullptr);
		}
	}
}

/**
 * @brief tests Bulge loop prediction
 * 
 * @date 2015-12-01
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
	int f =  0;
	for(int i = 0 ; i < 1000; i++)
	{
		// Initialize variables
		Sequence sequence = Sequence("GGGAAAGGGAAACCCCCC");
		std::string true_structure = "(((....((....)))))";
		
		Settings settings = Settings(0, nullptr, sequence);
		ReadData thermodynamics = ReadData();
		
		// Predict structure
		Zuker zuker = Zuker(settings, sequence, thermodynamics);
		zuker.energy();
		
		f++;
		test_pij(zuker);
		test_qij(zuker);
		test_vij(zuker);
		test_wij(zuker);
		test_nij2(zuker, sequence.size());
		test_loopmatrix(zuker);
		test_pathmatrix_corrected_from(zuker);
		
		
		zuker.traceback();
		
		// Obtain and compare results
		std::string predicted_structure;
		zuker.dot_bracket.format(sequence.size() , predicted_structure);
		
		bool valid_structure = predicted_structure.compare(true_structure) == 0;
		BOOST_CHECK_MESSAGE(valid_structure, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
		
		// The following code can be used to re-generate the funnctions above:
		/*
		test_vij
		for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
		{
			for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
			{
				Pair p = Pair(x,y);
				printf("	p = {%i, %i};e = zuker.vij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, %f),\"V(%i,%i) = \" + std::to_string(e) );\n", x,  y, zuker.vij.get(p), x, y);
			}
			printf("\n");
		}
		exit(1);
		
		test_wij
		/*for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
		{
			for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
			{
				Pair p = Pair(x,y);
				printf("	p = {%i, %i};e = zuker.wij.get(p);BOOST_CHECK_MESSAGE(IdenticalFloats(e, %f),\"W(%i,%i) = \" + std::to_string(e) );\n", x,  y, zuker.wij.get(p), x, y);
			}
			printf("\n");
		}
		exit(1);
		
		//test_loopmatrix
		for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
		{
			for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
			{
				Pair p = Pair(x,y);
				Pair p2 = zuker.loopmatrix.get(p);
				printf("	p = {%i, %i};p2 = zuker.loopmatrix.get(p);BOOST_CHECK_EQUAL(p2.first, %i);BOOST_CHECK_EQUAL(p2.second, %i);\n", x,  y, p2.first, p2.second);
			}
			printf("\n");
		}
		exit(1);
		
		//pathmatrix_corrected_from
		for(unsigned int x = 0;  x < (unsigned int) sequence.size(); x++)
		{
			for(unsigned int y = x;  y < (unsigned int) sequence.size(); y++)
			{
				Pair p = Pair(x,y);
				if( zuker.pathmatrix_corrected_from.get(p) )
				{
					printf("	p = {%i, %i};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == true ,\"path_matrix_corrected_from(%i,%i) = false\" );\n", x,  y, x, y);
				}
				else
				{
					printf("	p = {%i, %i};BOOST_CHECK_MESSAGE( zuker.pathmatrix_corrected_from.get(p) == false ,\"path_matrix_corrected_from(%i,%i) = true\" );\n", x,  y, x, y);
				}
			}
			printf("\n");
		}
		exit(1);
		*/
	}
}


BOOST_AUTO_TEST_SUITE_END()
