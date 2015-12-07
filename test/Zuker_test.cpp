/**
 * @file test/Zuker_test.cpp
 *
 * @date 2015-12-07
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



#define BOOST_TEST_MODULE Zuker


#include "main.hpp"

#include "Pair.hpp"
#include "Region.hpp"
#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"
#include "SubSequence.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"
#include "SegmentLoop.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"
#include "SegmentLoopTree.hpp"

#include "ScoringTree.hpp"

#include "ScoringMatrix.hpp"
#include "Settings.hpp"
#include "DotBracket.hpp"
#include "ReadData.hpp"

#include "Zuker.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Test_energy_loading) // Test energy loading


/**
 * @brief Tests predication of the canonical K-Turn (forward direction)
 *
 * @test
 *
 * @date 2015-08-06
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	Sequence rna_1 =   Sequence("guUGUGAUgaaacUGAac");
	size_t n = rna_1.size();
	
	std::string        segment_01_name  = "C/D-box K-turn";
	Sequence           segment_01_seq5p = Sequence("UGUGAU");
	std::vector <Pair> segment_01_bonds = { Pair({4, 1}), Pair({1, 1}), Pair({1, 1}) };
	Sequence           segment_01_seq3p = Sequence("UGA");// dbl check for reverse
	float              segment_01_nrg   = -100.0;
	Segment            segment_01       = Segment(segment_01_name, segment_01_seq5p, segment_01_bonds, segment_01_seq3p, segment_01_nrg);
	
	
	// sequence
	unsigned int i = 8;
	unsigned int j = 12;
	
	Nucleotide n1 = rna_1[i];
	Nucleotide n2 = rna_1[j];
	
	Pairing pairing = Pairing(n1, n2);
	
	
	// settings
	Sequence sequence = Sequence();
	ReadData thermodynamics = ReadData();
	
	//hairpin
	thermodynamics.tstackh[pairing.type][rna_1[i + 1]][rna_1[j - 1]] = -21.0;
	thermodynamics.triloop_map.clear();
	thermodynamics.triloop_map[Sequence("gaaac")] = -21.0;
	thermodynamics.loop_hairpin[3] = 0;
	
	// stem
	thermodynamics.stack[PairingType::GC][PairingType::UA] = -30.0;
	thermodynamics.stack[PairingType::UA][PairingType::UA] = -40.0;
	
	//k-turn
	thermodynamics.segments.insert(segment_01);
	
	
	// Init & run algorithm
	char *argv[] = { (char *) PACKAGE_NAME, (char *) "-s", (char *) "guUGUGAUgaaacUGAac" , NULL};
	int argc = sizeof(argv) / sizeof(char *) - 1;
	Settings settings = Settings(argc, argv, sequence);
	Zuker zuker = Zuker(settings, sequence , thermodynamics);
	
	BOOST_REQUIRE_EQUAL(zuker.energy() ,
						(float)
						- 100.0 +
						-21.0 +
						-21.0 +
						-30.0 +
						-1.50
					   );
					   
	zuker.traceback();
	
	std::string dotbracket_predicted;
	std::string dotbracket_valid = "((...((((...))))))";
	
	zuker.dot_bracket.format(n, dotbracket_predicted);
	
	BOOST_REQUIRE_EQUAL(dotbracket_predicted.compare(dotbracket_valid) , 0);
}



/**
 * @brief Tests prediction of the segmentloop
 *
 * @test
 *
 * @date 2015-12-07
 */
BOOST_AUTO_TEST_CASE(Test2_segmentloop)
{
	//5') gAA
	//    |||u
	//3') cAA
	Sequence rna_1 =   Sequence("gAAuAAc");
	size_t n = rna_1.size();
	
	std::string        segmentloop_01_name     = "k-loop test";
	Sequence           segmentloop_01_sequence = Sequence("AAuAA");
	std::vector <Pair> segmentloop_01_bonds    = { Pair({1, 1}), Pair({1, 1}) };
	float              segmentloop_01_nrg      = -100.0;
	SegmentLoop        segmentloop_01          = SegmentLoop(segmentloop_01_name, segmentloop_01_sequence , segmentloop_01_bonds, segmentloop_01_nrg);
	
	
	// subsequence
	unsigned int i = 1;
	unsigned int j = 5;
	
	Nucleotide n1 = rna_1[i];
	Nucleotide n2 = rna_1[j];
	
	Pairing pairing = Pairing(n1, n2);
	
	
	// settings
	Sequence sequence = Sequence();
	ReadData thermodynamics = ReadData();
	thermodynamics.loop_hairpin[5] = -5.0;
	
	//k-loop
	thermodynamics.segmentloops.insert(segmentloop_01);
	
	// Init & run algorithm
	char *argv[] = { (char *) PACKAGE_NAME, (char *) "-s", (char *) "cAAuAAg" , NULL};
	int argc = sizeof(argv) / sizeof(char *) - 1;
	Settings settings = Settings(argc, argv, sequence);
	Zuker zuker = Zuker(settings, sequence , thermodynamics);
	
	BOOST_CHECK_CLOSE_FRACTION(zuker.energy() , -100.0 + -2.10, 0.0001);
	
	zuker.traceback();
	
	std::string dotbracket_predicted;
	std::string dotbracket_valid = "(((.)))";
	
	zuker.dot_bracket.format(n, dotbracket_predicted);
	
	BOOST_REQUIRE_EQUAL(dotbracket_predicted.compare(dotbracket_valid) , 0);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(Test_matrices)


/**
 * @brief Check whether W contains the correct content for a not folding sequence
 *
 * @test Zuker::energy()
 *
 * @date 2015-07-15
 */
BOOST_AUTO_TEST_CASE(Test_W_matrix_01)
{
	// settings
	Sequence sequence = Sequence("aaa");
	ReadData thermodynamics = ReadData();
	
	// init
	char *argv[] = { (char *) PACKAGE_NAME, (char *) "-s", (char *) "aaa" , NULL};
	int argc = sizeof(argv) / sizeof(char *) - 1;
	
	Settings settings = Settings(argc, argv, sequence);
	
	// Run algorithm
	Zuker zuker = Zuker(settings, sequence , thermodynamics);
	
	unsigned int strlen = 3;
	unsigned int i, j;
	
	for(i = 0; i < strlen; i++)
	{
		for(j = i ; j < strlen; j++)
		{
			Pair p1 = Pair(i, j);
			
			BOOST_CHECK_EQUAL(zuker.wij.get(p1), 0.0);
			BOOST_CHECK_EQUAL(zuker.vij.get(p1), (i == j) ? N_INFINITY : 0.0);
		}
	}
	
	zuker.energy();
}


/**
 * @brief Tests if the non energy matrices are correctly filled
 *
 * @test
 *
 * @section DESCRIPTION
 * Prediction the following DotBracket structure: "(((...)))"
 *
 * The corresponding loopMatrix should be:
 *  <PRE>
 *  0,0 0,0 0,0 0,0 0,0 0,0 1,1 2,6 1,7
 *   .  0,0 0,0 0,0 0,0 0,0 2,2 2,6 2,6
 *   .   .  0,0 0,0 0,0 0,0 3,5 3,3 3,3
 *   .   .   .  0,0 0,0 0,0 0,0 0,0 0,0
 *   .   .   .   .  0,0 0,0 0,0 0,0 0,0
 *   .   .   .   .   .  0,0 0,0 0,0 0,0
 *   .   .   .   .   .   .  0,0 0,0 0,0
 *   .   .   .   .   .   .   .  0,0 0,0
 *   .   .   .   .   .   .   .   .  0,0
 *
 *  Pij must be:
 *  -2 -2 -2 -2  0  0  0  0 -1
 *   . -2 -2 -2 -2  1  1 -1  7
 *   .  . -2 -2 -2 -2 -1  6  6
 *   .  .  . -2 -2 -2 -2  3  3
 *   .  .  .  . -2 -2 -2 -2  4
 *   .  .  .  .  . -2 -2 -2 -2
 *   .  .  .  .  .  . -2 -2 -2
 *   .  .  .  .  .  .  . -2 -2
 *   .  .  .  .  .  .  .  . -2
 * </PRE>
 *
 * @date 2015-06-22
 *
 * @todo why do we see 1,1 and 2,2 while the other direction has 3,3?
 * @todo you can alsu use the scoring_matrixs' get_position(i,j) function to only store plain integers for memory efficiency
 */
BOOST_AUTO_TEST_CASE(Test_Sequence_GGGAAACCC)
{
	Sequence sequence = Sequence("GGGaaaCCC");
	
	char *argv[] = { (char *) PACKAGE_NAME, (char *) "-s", (char *) "GGGaaaCCC" , NULL};
	int argc = sizeof(argv) / sizeof(char *) - 1;
	
	unsigned int i = 2;
	unsigned int j = 6;
	
	Nucleotide p1 = sequence[i];
	Nucleotide p2 = sequence[j];
	
	Pairing pairing = Pairing(p1, p2);
	
	ReadData thermodynamics = ReadData();
	thermodynamics.tstackh[pairing.type][sequence[i + 1]][sequence[j - 1]] = -5.0;
	thermodynamics.triloop_map.clear();
	thermodynamics.triloop_map[Sequence("gaaac")] = -10.00;
	thermodynamics.loop_hairpin[3] = -7.5;
	thermodynamics.stack[PairingType::GC][PairingType::GC] = -5.0;
	
	Settings settings = Settings(argc, argv, sequence);
	
	// Run algorithm
	Zuker zuker = Zuker(settings, sequence , thermodynamics);
	float energy = zuker.energy();
	
	// loopmatrix
	// Checking these guys:                *
	//  -   -   -   -   -   -   -   -  1,7
	//      -   -   -   -   -   -  2,6  -
	//          -   -   -   -  3,5  -   -
	//              -   -  0,0  -   -   -
	//                  -   -   -   -   -
	//                      -   -   -   -
	//                          -   -   -
	//                              -   -
	//                                  -
	Pair pend = Pair(0, 8);
	Pair jump_1 = zuker.loopmatrix.get(pend);
	Pair jump_2 = zuker.loopmatrix.get(jump_1);
	Pair jump_3 = zuker.loopmatrix.get(jump_2);
	Pair jump_4 = zuker.loopmatrix.get(jump_3);
	
	BOOST_CHECK_EQUAL(jump_1.first == 1 , jump_1.second == 7);
	BOOST_CHECK_EQUAL(jump_2.first == 2 , jump_2.second == 6);
	BOOST_CHECK_EQUAL(jump_3.first == 3 , jump_3.second == 5);
	BOOST_CHECK_EQUAL(jump_4.first == 0 , jump_4.second == 0);
	
	// Path matrix
	// Checking these guys:        *
	//  -2 -2 -2 -2  0  0  0  0 -1
	//     -2 -2 -2 -2  1  1 -1  7
	//        -2 -2 -2 -2 -1  6  6
	//           -2 -2 -2 -2  3  3
	//              -2 -2 -2 -2  4
	//                 -2 -2 -2 -2
	//                    -2 -2 -2
	//                       -2 -2
	//                          -2
	BOOST_CHECK_EQUAL(zuker.pij.get(pend) , BOUND);
	BOOST_CHECK_EQUAL(zuker.pij.get(jump_1) , BOUND);
	BOOST_CHECK_EQUAL(zuker.pij.get(jump_2) , BOUND);
	BOOST_CHECK_EQUAL(zuker.pij.get(jump_3) , UNBOUND);
	
	// check total energy
	BOOST_CHECK_EQUAL(energy , -32.50);
	
	/*
	#if DEBUG
		zuker._print_loopmatrix(9);
		zuker._print_pij(9);
	#endif //DEBUG
	*/
}



/**
 * @brief Tests whether V and W succesfully sum up energy values
 *
 * @test Zuker::energy
 * @test Zuker::v
 * @test Zuker::wij
 *
 * @date 2015-07-13
 */
BOOST_AUTO_TEST_CASE(Test_vij_wij_01)
{
	Sequence sequence = Sequence("CGaaaaaCG");
	
	char *argv[] = { (char *) PACKAGE_NAME, (char *) "-s", (char *) "CGaaaaaCG" , NULL};
	int argc = sizeof(argv) / sizeof(char *) - 1;
	
	unsigned int i = 1;
	unsigned int j = 7;
	
	// Hairpin positions
	Nucleotide ph1 = sequence[i];
	Nucleotide ph2 = sequence[j];
	
	// Stack positions
	Nucleotide ps1 = sequence[i - 1];
	Nucleotide ps2 = sequence[j + 1];
	
	Pairing pairing_h = Pairing(ph1, ph2);
	Pairing pairing_s = Pairing(ps1, ps2);
	
	ReadData thermodynamics = ReadData();
	thermodynamics.tstackh[pairing_h.type][Nucleotide::A][Nucleotide::A] = -2.5;
	thermodynamics.triloop_map.clear();
	thermodynamics.loop_hairpin[5] = -2.5;
	thermodynamics.stack[PairingType::CG][PairingType::GC] = -10.0;
	
	// Run algorithm
	Settings settings = Settings(argc, argv, sequence);
	Zuker zuker = Zuker(settings, sequence , thermodynamics);
	
	// Give the hairpin loop a dGe of -100 (setting it explicitly in Vij), and set pij such that the position is considered to be calculated
	Pair p = Pair(i, j);
	zuker.pij.set(p, BOUND);
	zuker.vij.set(p, -100.0);
	
	/*
	#if DEBUG
		zuker._print_vij(9);
		zuker._print_wij(9);
	#endif //DEBUG
	*/
	
	float energy = zuker.energy();
	
	/*
	#if DEBUG
		zuker._print_vij(9);
		zuker._print_wij(9);
	#endif //DEBUG
	*/
	
	// Ensure the total dge is -100 (hairpin) + -10 (for the stack)
	Pair q = Pair(0, 8);
	
	BOOST_CHECK_EQUAL(zuker.vij.get(p), -100.0);
	BOOST_CHECK_EQUAL(zuker.vij.get(q), -110.0);
	
	/*
	@note not sure if this is due to setting vij or whether this is normal behavior
	BOOST_CHECK_EQUAL(zuker.wij.get(p), -0.0);
	*/
	BOOST_CHECK_EQUAL(zuker.wij.get(q), -110.0);
}



/**
 * @brief Tests whether bifurcation prediction works correctly
 *
 * @test Zuker::energy_bifurcation
 *
 * @section DESCRIPTION
 * vij:
 * <PRE>
 *   inf   inf   inf   inf   0.0   0.0 -22.5 -12.5 -12.5   0.0   0.0   0.0   0.0   0.0
 * -----   inf   inf   inf   inf -12.5   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
 * ----- -----   inf   inf   inf   inf   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
 * ----- ----- -----   inf   inf   inf   inf   0.0   0.0   0.0   0.0   0.0   0.0   0.0
 * ----- ----- ----- -----   inf   inf   inf   inf   0.0   0.0   0.0   0.0   0.0   0.0
 * ----- ----- ----- ----- -----   inf   inf   inf   inf   0.0   0.0   0.0   0.0 -12.5
 * ----- ----- ----- ----- ----- -----   inf   inf   inf   inf   0.0   0.0   0.0 -12.5
 * ----- ----- ----- ----- ----- ----- -----   inf   inf   inf   inf   0.0   0.0 -22.5
 * ----- ----- ----- ----- ----- ----- ----- -----   inf   inf   inf   inf -12.5   0.0
 * ----- ----- ----- ----- ----- ----- ----- ----- -----   inf   inf   inf   inf   0.0
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   inf   inf   inf   inf
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   inf   inf   inf
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   inf   inf
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   inf
 * </PRE>
 *
 * wij:
 * <PRE>
 *   0.0   0.0   0.0   0.0   0.0 -12.5 -22.5 -22.5 -22.5 -22.5 -22.5 -22.5 -35.0 -45.0
 * -----   0.0   0.0   0.0   0.0 -12.5 -12.5 -12.5 -12.5 -12.5 -12.5 -12.5 -25.0 -35.0
 * ----- -----   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0 -12.5 -22.5
 * ----- ----- -----   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0 -12.5 -22.5
 * ----- ----- ----- -----   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0 -12.5 -22.5
 * ----- ----- ----- ----- -----   0.0   0.0   0.0   0.0   0.0   0.0   0.0 -12.5 -22.5
 * ----- ----- ----- ----- ----- -----   0.0   0.0   0.0   0.0   0.0   0.0 -12.5 -22.5
 * ----- ----- ----- ----- ----- ----- -----   0.0   0.0   0.0   0.0   0.0 -12.5 -22.5
 * ----- ----- ----- ----- ----- ----- ----- -----   0.0   0.0   0.0   0.0 -12.5 -12.5
 * ----- ----- ----- ----- ----- ----- ----- ----- -----   0.0   0.0   0.0   0.0   0.0
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   0.0   0.0   0.0   0.0
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   0.0   0.0   0.0
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   0.0   0.0
 * ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----   0.0
 * </PRE>
 *
 * @date 2015-07-13
 */
BOOST_AUTO_TEST_CASE(Test_energy_bifuraction)
{
	Sequence sequence = Sequence("GGaaaCCCCaaaGG");
	
	// init
	char *argv[] = { (char *) PACKAGE_NAME, (char *) "-s", (char *) "GGaaaCCCCaaaGG" , NULL};
	int argc = sizeof(argv) / sizeof(char *) - 1;
	
	ReadData thermodynamics = ReadData();
	
	thermodynamics.tstackh[PairingType::GC][Nucleotide::A][Nucleotide::A] = -10;
	thermodynamics.tstackh[PairingType::CG][Nucleotide::A][Nucleotide::A] = -10;
	
	thermodynamics.stack[PairingType::CG][PairingType::CG] = -10.0;
	thermodynamics.stack[PairingType::GC][PairingType::GC] = -10.0;
	
	thermodynamics.triloop_map.clear();
	thermodynamics.loop_hairpin[3] = -2.5;
	
	// Run algorithm
	Settings settings = Settings(argc, argv, sequence);
	Zuker zuker = Zuker(settings, sequence , thermodynamics);
	float energy = zuker.energy();
	
	
	Pair pair01 = Pair(1, 5);
	Pair pair02 = Pair(8, 12);
	
	BOOST_CHECK_EQUAL(zuker.vij.get(pair01) , -12.5);
	BOOST_CHECK_EQUAL(zuker.wij.get(pair01) , -12.5);
	
	BOOST_CHECK_EQUAL(zuker.vij.get(pair02) , -12.5);
	BOOST_CHECK_EQUAL(zuker.wij.get(pair02) , -12.5);
	
	Pair pair03 = Pair(0, 6);
	Pair pair04 = Pair(7, 13);
	
	BOOST_CHECK_EQUAL(zuker.vij.get(pair03) , -22.5);
	BOOST_CHECK_EQUAL(zuker.wij.get(pair03) , -22.5);
	
	BOOST_CHECK_EQUAL(zuker.vij.get(pair04) , -22.5);
	BOOST_CHECK_EQUAL(zuker.wij.get(pair04) , -22.5);
	
	Pair pair05 = Pair(0, 13);
	
	BOOST_CHECK_EQUAL(zuker.vij.get(pair05) , -0.0);// not paired itself, but a bifurcation
	BOOST_CHECK_EQUAL(zuker.wij.get(pair05) , -45.0);// sum of branches (-22.5 * 2)
	BOOST_CHECK_EQUAL(energy , -45.0);
}

///@todo test function for Zuker::traceback, Zuker::traceback_pop and Zuker::traceback_push

BOOST_AUTO_TEST_SUITE_END()
