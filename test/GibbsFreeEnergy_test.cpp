/**
 * @file test/GibbsFreeEnergy_test.cpp
 *
 * @date 2016-01-21
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



#define BOOST_TEST_MODULE GibbsFreeEnergy



#include <array>
#include <vector>


#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Pair.hpp"
#include "Region.hpp"

#include "Segment.hpp"
#include "SegmentLoop.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"
#include "SegmentLoopTree.hpp"

#include "ReadData.hpp"
#include "ReadSegments.hpp"

#include "GibbsFreeEnergy.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing_hairpin_loop_element)

/**
 * @brief Tests hairpinloops with 4 unpaired nucleotides
 *
 * @test GibbsFreeEnergy::get_hairpin_loop_element
 * @test GibbsFreeEnergy::get_loop_hairpin
 * @test GibbsFreeEnergy::get_tloop
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestHairpinLoop1)
{
	Sequence rna_1 = Sequence("aaaGGGGACaaa");
	
	unsigned int i = 3;
	unsigned int j = 8;
	
	Pair p = Pair(i, j);
	
	Nucleotide p1 = rna_1[i];
	Nucleotide p2 = rna_1[j];
	
	Pairing pairing = Pairing(p1, p2);
	
	ReadData thermodynamics = ReadData();
	thermodynamics.tstackh[pairing.type][rna_1[i + 1]][rna_1[j - 1]] = -5;
	thermodynamics.tloop_map.clear();
	thermodynamics.tloop_map[Sequence("GGGGAC")] = -10;
	thermodynamics.loop_hairpin[4] = 0;
	
	GibbsFreeEnergy gfe_1 = GibbsFreeEnergy(rna_1, thermodynamics);
	
	BOOST_CHECK_EQUAL(gfe_1.get_hairpin_loop_element(p) , -15.0);
}



/**
 * @brief Tests hairpinloops with 3 unpaired nucleotides
 *
 * @test GibbsFreeEnergy::get_hairpin_loop_element
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestHairpinLoop2)
{
	Sequence rna_1 =  Sequence("aaaGUUACaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
	
	unsigned int i = 3;
	unsigned int j = 7;
	
	Pair p = Pair(i, j);
	
	Nucleotide n1 = rna_1[i];
	Nucleotide n2 = rna_1[j];
	
	Pairing pairing = Pairing(n1, n2);
	
	ReadData thermodynamics = ReadData();
	thermodynamics.tstackh[pairing.type][rna_1[i + 1]][rna_1[j - 1]] = -5;//thermodynamics.tstackh[rna_1[i]][rna_1[j]][rna_1[i + 1]][rna_1[j - 1]] = -5;
	thermodynamics.triloop_map.clear();
	thermodynamics.triloop_map[Sequence("GUUAC")] = -10;
	thermodynamics.loop_hairpin[3] = 0;
	
	GibbsFreeEnergy gfe_1 = GibbsFreeEnergy(rna_1, thermodynamics);
	
	BOOST_CHECK_EQUAL(gfe_1.get_hairpin_loop_element(p) , -15.0);
}



/**
 * @brief Tests hairpinloops with hairpinloop = 30 chars (fixed energy value)
 *
 * @test GibbsFreeEnergy::get_hairpin_loop_element
 * @test GibbsFreeEnergy::get_tstackh
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestHairpinLoop3)
{
	ReadData thermodynamics = ReadData();
	
	Sequence rna =  Sequence("CaaaaaaaaaAaaaaaaaaaAaaaaaaaaaAG");// 30 is perfect
	unsigned int i = 0;
	unsigned int j = 31;
	
	Pair p = Pair(i, j);
	
	Nucleotide n_i = rna[i];
	Nucleotide n_j = rna[j];
	
	Pairing pairing = Pairing(n_i, n_j);
	
	thermodynamics.tstackh[pairing.type][rna[i + 1]][rna[j - 1]] = 0.0;//thermodynamics.tstackh[rna[i]][rna[j]][rna[i + 1]][rna[j - 1]] = 0.0;
	
	thermodynamics.loop_hairpin[30] = -100.0;
	
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	BOOST_CHECK(gfe.get_hairpin_loop_element(p) == -100.0);
}



/**
 * @brief Tests hairpinloops with hairpinloop = 31 chars (energy = fixed)
 *
 * @test GibbsFreeEnergy::get_hairpin_loop_element
 *
 * @date 2016-01-21
 */
BOOST_AUTO_TEST_CASE(TestHairpinLoop4)
{
	Sequence rna = Sequence("CaaaaaaaaaAaaaaaaaaaAaaaaaaaaaAaG");// 30+1 should use interpolated extension
	
	unsigned int i = 0;
	unsigned int j = 32;
	
	Pair p = Pair(i, j);
	
	Nucleotide n_i = rna[i];
	Nucleotide n_j = rna[j];
	
	Pairing pairing = Pairing(n_i, n_j);
	
	ReadData thermodynamics = ReadData();
	
	thermodynamics.tstackh[pairing.type][rna[i + 1]][rna[j - 1]] = 0.0;
	
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	float interpolated_hairpin_loop = thermodynamics.loop_hairpin[30] + (thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float) 31 / 30));
	BOOST_CHECK_EQUAL(gfe.get_hairpin_loop_element(p) , interpolated_hairpin_loop);
}



/**
 * @brief Tests hairpinloops that are poly-C  & get_miscloop(MISCLOOP_C_HAIRPIN_OF_3) & get_miscloop(MISCLOOP_C_HAIRPIN_INTERCEPT) & get_miscloop(MISCLOOP_C_HAIRPIN_SLOPE)
 *
 * @test
 *
 * @date 2016-01-20
 */
BOOST_AUTO_TEST_CASE(TestHairpinLoop5)
{
	ReadData thermodynamics = ReadData();
	
	thermodynamics.miscloop[MISCLOOP_C_HAIRPIN_OF_3] = 3.0;
	thermodynamics.miscloop[MISCLOOP_C_HAIRPIN_INTERCEPT] = 4.0;
	thermodynamics.miscloop[MISCLOOP_C_HAIRPIN_SLOPE] = 4.0;
	
	Sequence rna_1 =  Sequence("AcccU");
	Sequence rna_2 =  Sequence("AccccU");
	Sequence rna_3 =  Sequence("AcccccU");
	
	unsigned int i = 0;
	unsigned int j_1 = 4;
	unsigned int j_2 = 5;
	unsigned int j_3 = 6;
	
	Pair p1 = Pair(i, j_1);
	Pair p2 = Pair(i, j_2);
	Pair p3 = Pair(i, j_3);
	
	Nucleotide n_i = rna_1[i];
	Nucleotide n_j1 = rna_1[j_1];
	Nucleotide n_j2 = rna_2[j_2];
	
	Pairing pairing1 = Pairing(n_i, n_j1);
	Pairing pairing2 = Pairing(n_i, n_j2);
	
	thermodynamics.tstackh[pairing1.type][rna_1[i + 1]][rna_1[j_1 - 1]] = 0.0;
	thermodynamics.tstackh[pairing2.type][rna_2[i + 1]][rna_2[j_2 - 1]] = 0.0;
	
	thermodynamics.loop_hairpin[3] = 0.0;
	thermodynamics.loop_hairpin[4] = 0.0;
	thermodynamics.loop_hairpin[5] = -2.0;
	
	GibbsFreeEnergy gfe_1 = GibbsFreeEnergy(rna_1, thermodynamics);
	GibbsFreeEnergy gfe_2 = GibbsFreeEnergy(rna_2, thermodynamics);
	GibbsFreeEnergy gfe_3 = GibbsFreeEnergy(rna_3, thermodynamics);
	
	BOOST_CHECK_EQUAL(gfe_1.get_hairpin_loop_element(p1), 3.0);
	BOOST_CHECK_EQUAL(gfe_2.get_hairpin_loop_element(p2), 20.0);
	BOOST_CHECK_EQUAL(gfe_3.get_hairpin_loop_element(p3), (6 * 4.0) - 2.0);
}



/**
 * @brief Tests hairpinloops have a GGG_U_loop & get_miscloop(MISCLOOP_GGG_U_PENALTY)
 *
 * @test GibbsFreeEnergy::get_hairpin_loop_element
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestHairpinLoop6)
{
	float energy_gggu_loop = -100.0;
	float energy_gu_pairing = 50;
	float energy_tloop_mask = 0.0;
	float energy_loop_hairpin = 0.0;
	
	// Test 1 data
	Sequence rna1 =  Sequence("cggGaaaaUccg");// sequence with GGG-U penalty
	
	unsigned int i1 = 3;
	unsigned int j1 = 8;
	
	Pair p1 = Pair(i1, j1);
	
	Nucleotide p1_i = rna1[i1];
	Nucleotide p1_j = rna1[j1];
	
	Pairing pairing1 = Pairing(p1_i, p1_j);
	
	// Test 2 data
	Sequence rna2 =  Sequence("ggGaaaaUcc");// sequence with GGG-U penalty
	
	unsigned int i2 = 2;
	unsigned int j2 = 7;
	
	Pair p2 = Pair(i2, j2);
	
	Nucleotide p2_i = rna2[i2];
	Nucleotide p2_j = rna2[j2];
	
	Pairing pairing2 = Pairing(p2_i, p2_j);
	
	// Test 3 data
	Sequence rna3 =  Sequence("gGaaaaUc");// sequence without GGG-U penalty
	
	unsigned int i3 = 1;
	unsigned int j3 = 6;
	
	Pair p3 = Pair(i3, j3);
	
	Nucleotide p3_i = rna3[i3];
	Nucleotide p3_j = rna3[j3];
	
	Pairing pairing3 = Pairing(p3_i, p3_j);
	
	// Parameter definitions
	ReadData thermodynamics = ReadData();
	
	thermodynamics.miscloop[MISCLOOP_GGG_U_PENALTY] = energy_gggu_loop;
	thermodynamics.tstackh[pairing1.type][rna1[i1 + 1]][rna1[j1 - 1]] = energy_gu_pairing;
	thermodynamics.tstackh[pairing2.type][rna2[i2 + 1]][rna2[j2 - 1]] = energy_gu_pairing;
	thermodynamics.tstackh[pairing3.type][rna3[i3 + 1]][rna3[j3 - 1]] = energy_gu_pairing;
	thermodynamics.tloop_map[Sequence("GaaaaU")] = energy_tloop_mask;
	thermodynamics.loop_hairpin[4] = energy_loop_hairpin;
	
	// Compare energy
	GibbsFreeEnergy gfe1 = GibbsFreeEnergy(rna1, thermodynamics);
	GibbsFreeEnergy gfe2 = GibbsFreeEnergy(rna2, thermodynamics);
	GibbsFreeEnergy gfe3 = GibbsFreeEnergy(rna3, thermodynamics);
	
	BOOST_CHECK_EQUAL(gfe1.get_hairpin_loop_element(p1) , (energy_gggu_loop  + energy_gu_pairing + energy_tloop_mask + energy_loop_hairpin));
	BOOST_CHECK_EQUAL(gfe2.get_hairpin_loop_element(p2) , (energy_gggu_loop  + energy_gu_pairing + energy_tloop_mask + energy_loop_hairpin));
	BOOST_CHECK_EQUAL(gfe3.get_hairpin_loop_element(p3) , (energy_gu_pairing + energy_tloop_mask + energy_loop_hairpin));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestingBulgeloop)


/**
 * @brief Tests the bulge loop size (n=4) penalty
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 * <PRE>
 *        aa
 *       a  a
 * 5') cgc  gca
 *     || \ || a
 * 3') gc  gcga
 * </PRE>
 *
 * Hairpin:
 * <PRE>
 * 5') ...Ca
 *          a
 * 3') ...Ga
 * </PRE>
 *
 * Bulge loop:
 * <PRE>
 * 5') ...CaaaaG...
 *         \  /
 * 3')   ...GC...
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_bulge_loop_element
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestBulgeLoop1)
{
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_CaaaG = 0.0;
	
	float energy_bulge_loopsize = -50.0;
	
	Sequence rna = Sequence("cgCaaaaGcaaagCGcg");
	
	// Hairpin
	unsigned int hi = 8;
	unsigned int hj = 12;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Bulge
	unsigned int i = 2;
	unsigned int j = 14;
	
	unsigned int ip = 7;
	unsigned int jp = 13;
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	
	ReadData thermodynamics = ReadData();
	
	// Ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("CaaaG")] = energy_hairpin_triloop_CaaaG;
	thermodynamics.loop_hairpin[3] = energy_loop_hairpin;
	
	// Bulge loop settings
	thermodynamics.loop_bulge[4] = energy_bulge_loopsize;
	
	// Compare energy
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	BOOST_CHECK_EQUAL(gfe.get_bulge_loop_element(r) ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_CaaaG +
					  energy_bulge_loopsize);
}



/**
 * @brief For all 1-length bulge loops, the Gibbs free energy is given with respect to the sequence of the surrounding stem.
 *
 * @section DESCRIPTION
 *     Predicts the Gibbs free energy of the following structure:
 *
 * <PRE>
 *           a
 *   5') cgc  gca
 *       || \ || a
 *   3') gc  gcga
 * </PRE>
 *
 * Hairpin:
 * <PRE>
 *  5') ...Ca
 *           a
 *  3') ...Ga
 * </PRE>
 *
 *  Bulge loop:
 * <PRE>
 *  5') ...CaG...
 *          \|
 *  3')  ...GC...
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_bulge_loop_element
 *
 * @date 2015-06-22
 */
BOOST_AUTO_TEST_CASE(TestBulgeLoop2)
{
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_CaaaG = 0.0;
	
	float energy_bulge_loopsize = -50.0;
	float energy_bulge_pairing_cggc = -25.0;
	
	Sequence rna = Sequence("cgCaGcaaagCGcg");
	
	// Hairpin
	unsigned int hi = 5;
	unsigned int hj = 9;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Bulge
	unsigned int i = 2;
	unsigned int j = 11;
	
	unsigned int ip = 4;
	unsigned int jp = 10;
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	ReadData thermodynamics = ReadData();
	
	// Ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("CaaaG")] = energy_hairpin_triloop_CaaaG;
	thermodynamics.loop_hairpin[3] = energy_loop_hairpin;
	
	// Bulge loop settings
	thermodynamics.stack[PairingType::CG][PairingType::GC] = energy_bulge_pairing_cggc;
	thermodynamics.loop_bulge[1] = energy_bulge_loopsize;
	
	// Compare energy
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	BOOST_CHECK_EQUAL(gfe.get_bulge_loop_element(r),
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_CaaaG +
					  energy_bulge_loopsize + energy_bulge_pairing_cggc);
}


/**
 * @brief Tests the bulge loop size penalty (32*a) and stack
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 * <PRE>
 *       aaaaaaaaaaaaaaaa
 *       a               a
 *       aa    aaaaaaaaaa
 *         a  a
 *   5') cgc  gca
 *       || \ || a
 *   3') gc  gcga
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_bulge_loop_element
 *
 * @date 2016-01-21
 */
BOOST_AUTO_TEST_CASE(TestBulgeLoop3)
{
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_CaaaG = 0.0;
	
	Sequence rna = Sequence("cgCAAAAAAAAAAaaaaaaaaaaAAAAAAAAAAaaGcaaagCGcg");
	
	// Hairpin
	unsigned int hi = 5 + 32 - 1;
	unsigned int hj = 9 + 32 - 1;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	Pairing hpairing = Pairing(hni, hnj);
	
	unsigned int i = 2;
	unsigned int j = 11 + 32 - 1;
	
	unsigned int ip = 4 + 32 - 1;
	unsigned int jp = 10 + 32 - 1;
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	
	ReadData thermodynamics = ReadData();
	
	// Ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("CaaaG")] = energy_hairpin_triloop_CaaaG;
	thermodynamics.loop_hairpin[3] = energy_loop_hairpin;
	
	// Compare energy
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);			// interpolation of bulge loop params should take place here.
	
	float intrapolated_energy = thermodynamics.loop_bulge[30] + (thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float) 32 / 30));
	
	BOOST_CHECK_EQUAL(gfe.get_bulge_loop_element(r) ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_CaaaG +
					  intrapolated_energy);
}



/**
 * @brief Tests the bulge loop size 4 and two * the AU penalty
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 * <PRE>
 *       aaaa
 * 5') cgA  Uca
 *     || \ || a
 * 3') gc  UAga
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_bulge_loop_element
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestBulgeLoop4)
{
	// Define energies
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_CaaaG = 0.0;
	
	float energy_AU_UA_penalty = -15.0;
	float energy_bulge_loopsize = -20.0;
	
	// Set sequence
	Sequence rna = Sequence("cgAaaaaUcaaagAUcg");
	
	// Hairpin
	unsigned int hairpin_size = 3;
	unsigned int hi = 8;
	unsigned int hj = hi + hairpin_size + 1;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::C);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::G);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Interior loop
	unsigned int i = 2;
	unsigned int j = 11 + 4 - 1;
	
	unsigned int ip = 4 + 4 - 1;
	unsigned int jp = 10 + 4 - 1;
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);			// interpolation of bulge loop params should take place here.
	
	// Set penalties: ensure the hairpin gets 0.0 Gibbs free energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("CaaaG")] = energy_hairpin_triloop_CaaaG;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: bulge loop penalties
	thermodynamics.miscloop[MISCLOOP_AU_PENALTY] = energy_AU_UA_penalty;
	thermodynamics.loop_bulge[4] = energy_bulge_loopsize;
	
	// Do check
	BOOST_CHECK_EQUAL(gfe.get_bulge_loop_element(r) ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_CaaaG +
					  energy_bulge_loopsize + (2 * energy_AU_UA_penalty));
}
///@todo final check: check whether the relative maximum bulge size fits with the length of the interpolated bulge loop penalty vector

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestingInteriorloop)

/**
 * @brief Tests assignment of Gibbs free energy for interior loop
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 * <PRE>
 *        aaa
 * 5') gcG   Cga
 *     |||   || a
 * 3') cgC   Gca
 *        aaa
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_interior_loop_element
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestInteriorLoop1)
{
	// Define energies
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_GaaaC = 0.0;
	float energy_interior_loopsize = -50.0;
	float energy_interior_GC_A_A = -10;
	
	// Set sequence
	Sequence rna = Sequence("gcGaaaCgaaacGaaaCgc");
	
	// Hairpin
	unsigned int hairpin_size = 3;
	unsigned int hi = 7;
	unsigned int hj = hi + hairpin_size + 1 ;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Interior loop
	unsigned int i = 2;
	unsigned int j = 16;
	
	unsigned int ip = 6;
	unsigned int jp = 12;
	
	Nucleotide ni = rna[i];
	Nucleotide nj = rna[j];
	
	Nucleotide nip = rna[ip];
	Nucleotide njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1], Nucleotide::A);
	
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);			// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.loop_interior[5] = N_INFINITY;
	thermodynamics.loop_interior[6] = energy_interior_loopsize;
	thermodynamics.loop_interior[7] = N_INFINITY;
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	
	// Check
	BOOST_CHECK_EQUAL(gfe.get_interior_loop_element(r) ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  energy_interior_loopsize + (2 * energy_interior_GC_A_A));
}



/**
 * @brief Tests energy assignment to Interior Loop
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 *
 * <PRE>
 *       {r1}
 * 5') gcG  Cga
 *     |||  || a
 * 3') cgC  Gca
 *       {r2}
 * </PRE>
 *
 * Where for:
 *
 * n_unpaired=34:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaAA
 * {r2} = AaaaaAaaaaAaaaaAA
 * </PRE>
 *
 * n_unpaired=35:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaAAA
 * {r2} = AaaaaAaaaaAaaaaAA
 * </PRE>
 *
 * n_unpaired=36:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaAAA
 * {r2} = AaaaaAaaaaAaaaaAAA
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_interior_loop_element
 *
 * @date 2016-01-21
 */
BOOST_AUTO_TEST_CASE(TestInteriorLoop2)
{
	// Define energies
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_GaaaC = 0.0;
	float energy_interior_GC_A_A = -10;
	
	float asymetric_loopsize_correction = 0;
	
	// ------------------------ n_unpaired = 34 ------------------------
	// Set sequence
	Sequence rna = Sequence("gcGAaaaaAaaaaAaaaaAACgaaacGAaaaaAaaaaAaaaaAACgc");//n_unpaired = AaaaaAaaaaAaaaaAA + AaaaaAaaaaAaaaaAA = 34
	
	// Hairpin
	unsigned int hairpin_size = 3;
	
	// Hairpin definition for interior loop of size n=34 block
	unsigned int hi = 21;
	unsigned int hj = hi + hairpin_size + 1 ;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Interior loop
	unsigned int n_unpaired_r1 = 17;
	unsigned int n_unpaired_r2 = 17;
	
	unsigned int i = 2;
	unsigned int j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	unsigned int ip = i + n_unpaired_r1 + 1;
	unsigned int jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	Nucleotide ni = rna[i];
	Nucleotide nj = rna[j];
	
	Nucleotide nip = rna[ip];
	Nucleotide njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1], Nucleotide::A);
	
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe_34 = GibbsFreeEnergy(rna, thermodynamics);		// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	thermodynamics.miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] = asymetric_loopsize_correction;
	
	// Check
	float energy_34 = gfe_34.get_interior_loop_element(r);
	float intrapolated_energy = thermodynamics.loop_interior[30] +
								(thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float)(n_unpaired_r1 + n_unpaired_r2) / 30));
								
	BOOST_CHECK_EQUAL(energy_34 ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  intrapolated_energy + (2 * energy_interior_GC_A_A));
					  
					  
	// ------------------------ n_unpaired = 35 ------------------------
	// Set sequence
	rna = Sequence("gcGAaaaaAaaaaAaaaaAAACgaaacGAaaaaAaaaaAaaaaAACgc");//n_unpaired = AaaaaAaaaaAaaaaAAA + AaaaaAaaaaAaaaaAA = 35
	
	// Hairpin definition for interior loop of size n=35 block
	hi = 22;
	hj = hi + hairpin_size + 1 ;
	
	hni = rna[hi];
	hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	hpairing = Pairing(hni, hnj);
	
	// Interior loop
	n_unpaired_r1 = 18;
	n_unpaired_r2 = 17;
	
	i = 2;
	j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	ip = i + n_unpaired_r1 + 1;
	jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	ni = rna[i];
	nj = rna[j];
	
	nip = rna[ip];
	njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1] , Nucleotide::A);
	
	
	p = Pair {i, j};
	pp = Pair {ip, jp};
	
	Region r2 = Region {p, pp};
	
	// Compare energy
	thermodynamics = ReadData();
	GibbsFreeEnergy gfe_35 = GibbsFreeEnergy(rna, thermodynamics);		// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	thermodynamics.miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] = asymetric_loopsize_correction;
	
	// Check
	float energy_35 = gfe_35.get_interior_loop_element(r2);
	intrapolated_energy = thermodynamics.loop_interior[30] +
						  (thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float)(n_unpaired_r1 + n_unpaired_r2) / 30));
						  
	BOOST_CHECK_EQUAL(energy_35 ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  intrapolated_energy + (2 * energy_interior_GC_A_A));
					  
					  
	// ------------------------ n_unpaired = 36 ------------------------
	// Set sequence
	rna = Sequence("gcGAaaaaAaaaaAaaaaAAACgaaacGAaaaaAaaaaAaaaaAAACgc");////n_unpaired = AaaaaAaaaaAaaaaAAA + AaaaaAaaaaAaaaaAAA = 36
	
	// Hairpin definition for interior loop of size n=36 block
	hi = 22;
	hj = hi + hairpin_size + 1 ;
	
	hni = rna[hi];
	hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1], Nucleotide::A);
	
	
	hpairing = Pairing(hni, hnj);
	
	// Interior loop
	n_unpaired_r1 = 18;
	n_unpaired_r2 = 18;
	
	i = 2;
	j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	ip = i + n_unpaired_r1 + 1;
	jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	ni = rna[i];
	nj = rna[j];
	
	nip = rna[ip];
	njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	p = Pair {i, j};
	pp = Pair {ip, jp};
	
	Region r3 = Region {p, pp};
	
	// Compare energy
	thermodynamics = ReadData();
	GibbsFreeEnergy gfe_36 = GibbsFreeEnergy(rna, thermodynamics);		// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	thermodynamics.miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] = asymetric_loopsize_correction;
	
	// Check
	float energy_36 = gfe_36.get_interior_loop_element(r3);
	intrapolated_energy = thermodynamics.loop_interior[30] +
						  (thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float)(n_unpaired_r1 + n_unpaired_r2) / 30));
						  
	BOOST_CHECK_EQUAL(energy_36 ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  intrapolated_energy + (2 * energy_interior_GC_A_A));
					  
	// -----------------------------------------------------------------
	
	BOOST_CHECK(energy_34 < energy_35);
	BOOST_CHECK(energy_35 < energy_36);
}



/**
 * @brief Checks wheter the poppen function and limitation using miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] work.
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 * <PRE>
 *       {r1}
 * 5') gcG  Cga
 *     |||  || a
 * 3') cgC  Gca
 *       {r2}
 * </PRE>
 *
 * The poppen vector contains only 4 values and therefore we will test it with four values as well.
 * Where {r1} and {r2} are stretches of sequences that remain unpaired, and for:
 *
 * n_unpaired=32:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaA
 * {r2} = AaaaaAaaaaAaaaaA
 * </PRE>
 *
 * n_unpaired=33:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaA
 * {r2} = AaaaaAaaaaAaaaaAA
 * </PRE>
 *
 * n_unpaired=34:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaA
 * {r2} = AaaaaAaaaaAaaaaAAA
 * </PRE>
 *
 * n_unpaired=35:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaA
 * {r2} = AaaaaAaaaaAaaaaAAAA
 * </PRE>
 *
 * n_unpaired=36:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaA
 * {r2} = AaaaaAaaaaAaaaaAAAAA
 * </PRE>
 *
 * n_unpaired=37:
 * <PRE>
 * {r1} = AaaaaAaaaaAaaaaA
 * {r2} = AaaaaAaaaaAaaaaAAAAAA
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_interior_loop_element
 *
 * @date 2016-01-21
 */
BOOST_AUTO_TEST_CASE(TestInteriorLoop2_poppen)
{
	// Define energies
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_GaaaC = 0.0;
	
	float energy_interior_GC_A_A = -10.0;
	
	float asymetric_loopsize_correction = 2.0;
	float asymetric_loopsize_poppen_4 = 1.0;
	
	// ------------------------ n_unpaired = 32 ------------------------
	// Set sequence
	Sequence rna = Sequence("gcGAaaaaAaaaaAaaaaACgaaacGAaaaaAaaaaAaaaaACgc");//n_unpaired = AaaaaAaaaaAaaaaA + AaaaaAaaaaAaaaaA = 32
	
	// Hairpin
	unsigned int hairpin_size = 3;
	
	// Hairpin definition for interior loop of size n=32 block
	unsigned int hi = 20;
	unsigned int hj = hi + hairpin_size + 1 ;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Interior loop
	unsigned int n_unpaired_r1 = 16;
	unsigned int n_unpaired_r2 = 16;
	
	unsigned int i = 2;
	unsigned int j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	unsigned int ip = i + n_unpaired_r1 + 1;
	unsigned int jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	Nucleotide ni = rna[i];
	Nucleotide nj = rna[j];
	
	Nucleotide nip = rna[ip];
	Nucleotide njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1], Nucleotide::A);
	
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe_32 = GibbsFreeEnergy(rna, thermodynamics);		// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	thermodynamics.miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] = asymetric_loopsize_correction;
	thermodynamics.poppen_p[1] = 0.1f;
	thermodynamics.poppen_p[2] = 0.2f;
	thermodynamics.poppen_p[3] = 0.3f;
	thermodynamics.poppen_p[4] = asymetric_loopsize_poppen_4;
	
	// Check
	float energy_32 = gfe_32.get_interior_loop_element(r);
	float intrapolated_energy = thermodynamics.loop_interior[30] +
								(thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float)(n_unpaired_r1 + n_unpaired_r2) / 30));
								
	BOOST_CHECK_EQUAL(energy_32 ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  intrapolated_energy + (2 * energy_interior_GC_A_A) +
					  (0 * asymetric_loopsize_poppen_4)
					 );
					 
					 
	// ------------------------ n_unpaired = 33 ------------------------
	// Set sequence
	rna = Sequence("gcGAaaaaAaaaaAaaaaACgaaacGAaaaaAaaaaAaaaaAACgc");//n_unpaired = AaaaaAaaaaAaaaaA + AaaaaAaaaaAaaaaAA = 33
	
	// Hairpin definition for interior loop of size n=33 block
	hi = 20;
	hj = hi + hairpin_size + 1 ;
	
	hni = rna[hi];
	hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	hpairing = Pairing(hni, hnj);
	
	// Interior loop
	n_unpaired_r1 = 16;
	n_unpaired_r2 = 17;
	
	i = 2;
	j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	ip = i + n_unpaired_r1 + 1;
	jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	ni = rna[i];
	nj = rna[j];
	
	nip = rna[ip];
	njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1] , Nucleotide::A);
	
	
	p = Pair {i, j};
	pp = Pair {ip, jp};
	
	Region r2 = Region {p, pp};
	
	// Compare energy
	GibbsFreeEnergy gfe_33 = GibbsFreeEnergy(rna, thermodynamics);		// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	thermodynamics.miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] = asymetric_loopsize_correction;
	thermodynamics.poppen_p[1] = 0.1f;
	thermodynamics.poppen_p[2] = 0.2f;
	thermodynamics.poppen_p[3] = 0.3f;
	thermodynamics.poppen_p[4] = asymetric_loopsize_poppen_4;
	
	// Check
	float energy_33 = gfe_33.get_interior_loop_element(r2);
	intrapolated_energy = thermodynamics.loop_interior[30] +
						  (thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float)(n_unpaired_r1 + n_unpaired_r2) / 30));
						  
	BOOST_CHECK_EQUAL(energy_33 ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  intrapolated_energy + (2 * energy_interior_GC_A_A) +
					  (1 * asymetric_loopsize_poppen_4)
					 );
	// ------------------------ n_unpaired = 34 ------------------------
	// Set sequence
	rna = Sequence("gcGAaaaaAaaaaAaaaaACgaaacGAaaaaAaaaaAaaaaAAACgc");//n_unpaired = AaaaaAaaaaAaaaaA + AaaaaAaaaaAaaaaAAA = 34
	
	// Hairpin definition for interior loop of size n=34 block
	hi = 20;
	hj = hi + hairpin_size + 1 ;
	
	hni = rna[hi];
	hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	hpairing = Pairing(hni, hnj);
	
	// Interior loop
	n_unpaired_r1 = 16;
	n_unpaired_r2 = 18;
	
	i = 2;
	j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	ip = i + n_unpaired_r1 + 1;
	jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	ni = rna[i];
	nj = rna[j];
	
	nip = rna[ip];
	njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1] , Nucleotide::A);
	
	
	p = Pair {i, j};
	pp = Pair {ip, jp};
	
	Region r3 = Region {p, pp};
	
	// Compare energy
	thermodynamics = ReadData();
	GibbsFreeEnergy gfe_34 = GibbsFreeEnergy(rna, thermodynamics);		// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	thermodynamics.miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] = asymetric_loopsize_correction;
	thermodynamics.poppen_p[1] = 0.1f;
	thermodynamics.poppen_p[2] = 0.2f;
	thermodynamics.poppen_p[3] = 0.3f;
	thermodynamics.poppen_p[4] = asymetric_loopsize_poppen_4;
	
	// Check
	float energy_34 = gfe_34.get_interior_loop_element(r3);
	intrapolated_energy = thermodynamics.loop_interior[30] +
						  (thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float)(n_unpaired_r1 + n_unpaired_r2) / 30));
						  
	BOOST_CHECK_EQUAL(energy_34 ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  intrapolated_energy + (2 * energy_interior_GC_A_A) +
					  (2 * asymetric_loopsize_poppen_4)
					 );
					 
					 
	// ------------------------ n_unpaired = 35 ------------------------
	// Set sequence
	rna = Sequence("gcGAaaaaAaaaaAaaaaACgaaacGAaaaaAaaaaAaaaaAAAACgc");//n_unpaired = AaaaaAaaaaAaaaaA + AaaaaAaaaaAaaaaAAAA = 35
	
	// Hairpin definition for interior loop of size n=35 block
	hi = 20;
	hj = hi + hairpin_size + 1 ;
	
	hni = rna[hi];
	hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	hpairing = Pairing(hni, hnj);
	
	// Interior loop
	n_unpaired_r1 = 16;
	n_unpaired_r2 = 19;
	
	i = 2;
	j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	ip = i + n_unpaired_r1 + 1;
	jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	ni = rna[i];
	nj = rna[j];
	
	nip = rna[ip];
	njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1] , Nucleotide::A);
	
	
	p = Pair {i, j};
	pp = Pair {ip, jp};
	
	Region r4 = Region {p, pp};
	
	// Compare energy
	thermodynamics = ReadData();
	GibbsFreeEnergy gfe_35 = GibbsFreeEnergy(rna, thermodynamics);		// interpolation of interior loop params (usually n_nunpaired > 30) should take place here.
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.tstki[PairingType::GC][Nucleotide::A][Nucleotide::A] = energy_interior_GC_A_A;
	thermodynamics.miscloop[MISCLOOP_ASYMETRIC_INTENRAL_LOOP] = asymetric_loopsize_correction;
	thermodynamics.poppen_p[1] = 0.1f;
	thermodynamics.poppen_p[2] = 0.2f;
	thermodynamics.poppen_p[3] = 0.3f;
	thermodynamics.poppen_p[4] = asymetric_loopsize_poppen_4;
	
	// Check
	float energy_35 = gfe_35.get_interior_loop_element(r4);
	intrapolated_energy = thermodynamics.loop_interior[30] +
						  (thermodynamics.miscloop[MISCLOOP_PRELOG] * (float) log((float)(n_unpaired_r1 + n_unpaired_r2) / 30));
						  
	BOOST_CHECK_EQUAL(energy_35 ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  intrapolated_energy + (2 * energy_interior_GC_A_A) +
					  asymetric_loopsize_correction
					 );
					 
	// -----------------------------------------------------------------
}



/**
 * @brief Tests the energy assigned to an interior loop element
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 * <PRE>
 * 5') [i] [i+1] [i'] ...
 *      |         |
 * 3') [j] [j-1] [j'] ...
 * </PRE>
 *
 * @test get_int11 & get_interior_loop_element
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestInteriorLoop3)
{
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_GaaaC = 0.0;
	float energy_interior_loop_penalty2 = -75.0;
	
	Sequence rna = Sequence("gcGaCgaaacGaCgc");
	
	// Check hairpin
	unsigned int hairpin_size = 3;
	
	unsigned int hi = 5;
	unsigned int hj = hi + hairpin_size + 1 ;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Check interior loop
	unsigned int n_unpaired_r1 = 1;
	unsigned int n_unpaired_r2 = 1;
	
	unsigned int i = 2;
	unsigned int j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	unsigned int ip = i + n_unpaired_r1 + 1;
	unsigned int jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	Nucleotide ni = rna[i];
	Nucleotide nj = rna[j];
	
	Nucleotide nip = rna[ip];
	Nucleotide njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1], Nucleotide::A);
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.int11[PairingType::GC][PairingType::CG][Nucleotide::A][Nucleotide::A] =  energy_interior_loop_penalty2;
	
	// Check
	float energy = gfe.get_interior_loop_element(r);
	
	BOOST_CHECK_EQUAL(energy ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  energy_interior_loop_penalty2
					 );
}



/**
 * @brief Tests get_int21() & get_interior_loop_element()
 *
 * @section DESCRIPTION
 * Predicts the Gibbs free energy of the following structure:
 * <PRE>
 * 5') [i] [i+1]       [i'] ...
 *      |               |
 * 3') [j] [j-1] [j-2] [j'] ...
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_interior_loop_element
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestInteriorLoop4)
{
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_GaaaC = 0.0;
	float energy_interior_loop_penalty2 = -75.0;
	
	Sequence rna = Sequence("gcGaCgaaacGaaCgc");
	
	// Check hairpin
	unsigned int hairpin_size = 3;
	
	unsigned int hi = 5;
	unsigned int hj = hi + hairpin_size + 1 ;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Check interior loop
	unsigned int n_unpaired_r1 = 1;
	unsigned int n_unpaired_r2 = 2;
	
	unsigned int i = 2;
	unsigned int j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	unsigned int ip = i + n_unpaired_r1 + 1;
	unsigned int jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	Nucleotide ni = rna[i];
	Nucleotide nj = rna[j];
	
	Nucleotide nip = rna[ip];
	Nucleotide njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 1], Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[jp + 2] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[j - 2] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[jp + 1], Nucleotide::A);
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.int21
	[PairingType::GC]
	[PairingType::CG]
	[Nucleotide::A]
	[Nucleotide::A]
	[Nucleotide::A] = energy_interior_loop_penalty2;
	
	// Check
	float energy = gfe.get_interior_loop_element(r);
	
	BOOST_CHECK_EQUAL(energy ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  energy_interior_loop_penalty2
					 );
}



/**
 * @brief Tests whether the correct enegy is found for interior loop
 *
 * @section DESCRIPTION
 * Tests whether the correct enegy is found for interior loop of the following structure:
 * <PRE>
 * 5') [i] [i+1] [i+2] [i'] ...
 *      |               |
 * 3') [j] [j-1]       [j'] ...
 * </PRE>
 *
 * @test GibbsFreeEnergy::get_interior_loop_element
 * @test GibbsFreeEnergy::get_int21 (using reverse sequence)
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestInteriorLoop5)
{
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_GaaaC = 0.0;
	float energy_interior_loop_penalty2 = -75.0;
	
	Sequence rna = Sequence("gcGagCgaaacGaCgc");
	
	// Check hairpin
	unsigned int hairpin_size = 3;
	
	unsigned int hi = 6;
	unsigned int hj = hi + hairpin_size + 1 ;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Check interior loop
	unsigned int n_unpaired_r1 = 2;
	unsigned int n_unpaired_r2 = 1;
	
	unsigned int i = 2;
	unsigned int j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	unsigned int ip = i + n_unpaired_r1 + 1;
	unsigned int jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	Nucleotide ni = rna[i];
	Nucleotide nj = rna[j];
	
	Nucleotide nip = rna[ip];
	Nucleotide njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[ip - 2], Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[jp + 1] , Nucleotide::A);
	
	BOOST_CHECK_EQUAL(rna[i + 2] , Nucleotide::G);
	BOOST_CHECK_EQUAL(rna[ip - 1], Nucleotide::G);
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.int21
	[PairingType::GC]
	[PairingType::CG]
	[Nucleotide::A]
	[Nucleotide::G]
	[Nucleotide::A] =  energy_interior_loop_penalty2;
	
	// Check
	float energy = gfe.get_interior_loop_element(r);
	
	BOOST_CHECK_EQUAL(energy ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  energy_interior_loop_penalty2
					 );
}



/**
 * @brief Tests get_int22() & get_interior_loop_element()
 *
 * @section DESCRIPTION
 *    Predicts the Gibbs free energy of the following structure:
 *
 * <PRE>
 *   5') [i] [i+1] [i+2] [i'] ...
 *        |               |
 *   3') [j] [j-1] [j-2] [j'] ...
 * </PRE>
 *
 *  With interior loop sequence:
 * <PRE>
 *        aa
 * 5') gcG  Cga
 *     |||  || a
 * 3') cgC  Gca
 *        aa
 * </PRE>
 *
 * @test
 *
 * @date 2015-06-23
 */
BOOST_AUTO_TEST_CASE(TestInteriorLoop6)
{
	float energy_loop_hairpin = 0.0;
	float energy_hairpin_stack = 0.0;
	float energy_hairpin_triloop_GaaaC = 0.0;
	float energy_interior_loop_penalty2 = -75.0;
	
	Sequence rna = Sequence("gcGaaCgaaacGaaCgc");
	
	// Check hairpin
	unsigned int hairpin_size = 3;
	
	unsigned int hi = 6;
	unsigned int hj = hi + hairpin_size + 1 ;
	
	Nucleotide hni = rna[hi];
	Nucleotide hnj = rna[hj];
	
	BOOST_CHECK_EQUAL(hni , Nucleotide::G);
	BOOST_CHECK_EQUAL(hnj , Nucleotide::C);
	
	Pairing hpairing = Pairing(hni, hnj);
	
	// Check interior loop
	unsigned int n_unpaired_r1 = 2;
	unsigned int n_unpaired_r2 = 2;
	
	unsigned int i = 2;
	unsigned int j = i + n_unpaired_r1 + 2 + hairpin_size + 2 + n_unpaired_r2 + 1;
	
	unsigned int ip = i + n_unpaired_r1 + 1;
	unsigned int jp = i + n_unpaired_r1 + 2 + hairpin_size + 2;
	
	Nucleotide ni = rna[i];
	Nucleotide nj = rna[j];
	
	Nucleotide nip = rna[ip];
	Nucleotide njp = rna[jp];
	
	BOOST_CHECK_EQUAL(ni , Nucleotide::G);
	BOOST_CHECK_EQUAL(nj , Nucleotide::C);
	
	BOOST_CHECK_EQUAL(nip , Nucleotide::C);
	BOOST_CHECK_EQUAL(njp , Nucleotide::G);
	
	
	BOOST_CHECK_EQUAL(rna[i + 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 1] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[i + 2] , Nucleotide::A);
	BOOST_CHECK_EQUAL(rna[j - 2] , Nucleotide::A);
	
	
	Pair p = Pair {i, j};
	Pair pp = Pair {ip, jp};
	
	Region r = Region {p, pp};
	
	// Compare energy
	ReadData thermodynamics = ReadData();
	GibbsFreeEnergy gfe = GibbsFreeEnergy(rna, thermodynamics);
	
	// Set penalties: ensure the hairpin gets 0 energy
	thermodynamics.tstackh[hpairing.type][rna[hi + 1]][rna[hj - 1]] = energy_hairpin_stack;
	thermodynamics.tloop_map[Sequence("GaaaC")] = energy_hairpin_triloop_GaaaC;
	thermodynamics.loop_hairpin[hairpin_size] = energy_loop_hairpin;
	
	// Set penalties: interior loop penalties
	thermodynamics.int22
	[PairingType::GC]
	[PairingType::CG]
	[Nucleotide::A]
	[Nucleotide::A]
	[Nucleotide::A]
	[Nucleotide::A] = energy_interior_loop_penalty2;
	
	// Check
	float energy = gfe.get_interior_loop_element(r);
	
	BOOST_CHECK_EQUAL(energy ,
					  energy_loop_hairpin + energy_hairpin_stack + energy_hairpin_triloop_GaaaC +
					  energy_interior_loop_penalty2
					 );
}


///@todo final check: check whether the relative maximum hairpin size fits with the length of the interpolated hairpin loop penalty vector

///@todo test GibbsFreeEnergy::get_stacking_pair_element(Pair &arg_pair)

///@todo test Check hairpin size - at the moment it's stored in the settings AS WELL as in readdata << one of them should be removed

BOOST_AUTO_TEST_SUITE_END()
