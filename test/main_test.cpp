/**
 * @file test/main_test.cpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether an unstructured RNA is indeed unstructured
 *
 * @test
 *
 */
BOOST_AUTO_TEST_CASE(Test_unfolded)
{
	// Initialize variables
	Sequence sequence = Sequence("AAAaaaAAA");
	std::string true_structure = ".........";
	
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure); ///@todo make it size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
}



/**
 * @brief tests Hairpin sequence GGGAAACCC to be folded as (((...)))
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_hairpin)
{
	// Initialize variables
	Sequence sequence = Sequence("GGGaaaCCC");
	std::string true_structure = "(((...)))";
	
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	thermodynamics.triloop_map[Sequence("gaaac")] = -10.00;
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure); ///@todo make it size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
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
 *
 */
BOOST_AUTO_TEST_CASE(Test_bulge_loop)
{
	// Initialize variables
	Sequence sequence = Sequence("GGGAAAGGGAAACCCCCC");
	std::string true_structure = "(((...(((...))))))";
	
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	thermodynamics.triloop_map[Sequence("gaaac")] = -10.00;
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure); ///@size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
}

/**
 * @brief tests Interior loop prediction
 *
 * @test
 *
 * @section DESCRIPTION
 * Predicts folding of the following interior-loop structure:
 * <PRE>
 * 5') GGGGAAAGGGGA
 *     ||||   |||| A
 * 3') CCCCAAACCCCA
 * </PRE>
 *
 *  With the corresponding Sequence and DotBracket structure:
 * <PRE>
 * GGGGAAAGGGGAAACCCCAAACCCC
 * ((((...((((...))))...))))
 * </PRE>
 *
 */
BOOST_AUTO_TEST_CASE(Test_interior_loop)
{
	// Initialize variables
	Sequence sequence = Sequence("GGGGAAAGGGGAAACCCCAAACCCC");
	std::string true_structure = "((((...((((...))))...))))";
	
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	thermodynamics.triloop_map[Sequence("gaaac")] = -10.00;
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure); ///@todo unsinged int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
}

/**
 * @brief tests Bifurcation prediction
 *
 * @test
 *
 * @section DESCRIPTION
 * Predicts folding of the following bifurcation:
 * <PRE>
 *        A A
 *       G. A
 *      G. C
 *     G. C
 * GGGG  C
 * ||||
 * CCCC .G
 *     C .G
 *      C .G
 *       C  A
 *        A A
 * </PRE>
 *
 * With the corresponding Sequence and DotBracket structure:
 * <PRE>
 * GGGGgggaaaCCCgggAAAcccCCCC
 * (((((((...)))(((...)))))))
 * </PRE>
 *
 */
BOOST_AUTO_TEST_CASE(Test_bifurcation)
{
	// Initialize variables
	Sequence sequence = Sequence("GGGGgggaaaCCCgggAAAcccCCCC");
	std::string true_structure = "(((((((...)))(((...)))))))";
	
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	thermodynamics.triloop_map[Sequence("gaaac")] = -10.00;
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure); ///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
}



/**
 * @brief Runs the function al tests of all the segments in the test-file
 *
 * @test
 *
 * @section DESCRIPTION
 * Uses the same segment_file as it would use normally. In case this
 * functional test is executed with a modified version of this file
 * this might result into unneccesairy errors.
 *
 * @todo BOOST_REQUIRE_EQUAL << md5sum , segment_file
 */
BOOST_AUTO_TEST_CASE(Test_kturns)
{
	Sequence dummy = Sequence("A");
	Settings settings = Settings(0, nullptr, dummy);
	
	ReadData thermodynamics = ReadData();
	std::vector<rna_example> rna_examples;
	ReadSegments r = ReadSegments(settings.segment_filename);
	r.parse(thermodynamics.segments, thermodynamics.segmentloops, rna_examples);
	
	DotBracket db = DotBracket();
	
	// Test each example separately
	for(std::vector<rna_example>::iterator example = rna_examples.begin(); example != rna_examples.end(); ++example)
	{
		Zuker zuker = Zuker(settings, (*example).sequence, thermodynamics);
		zuker.energy();
		zuker.traceback();
		
		std::string predicted_structure = "";
		zuker.dot_bracket.format((unsigned int)(*example).sequence.size() , predicted_structure);  ///@todo unsigned int -> size_t
		
		BOOST_REQUIRE_EQUAL((*example).dot_bracket_pattern.size() , predicted_structure.size());
		BOOST_CHECK_MESSAGE(db.match((*example).dot_bracket_pattern, predicted_structure), "Predicted structure of '" << (*example).title << "' doesn't match it's true structure:\n\t" << predicted_structure << " (predicted structure)\n\t" << (*example).dot_bracket_pattern << " (pattern of true structure)\n");
	}
}



/**
 * @brief Runs the function al tests of all the segments in the test-file with segmentation-funcitonality disabled.
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_kturns_segments_disabled)
{
	Sequence dummy = Sequence("A");
	Settings settings = Settings(0, nullptr, dummy);
	settings.segment_prediction_functionality = false;
	
	ReadData thermodynamics = ReadData();
	std::vector<rna_example> rna_examples;
	
	DotBracket db = DotBracket();
	
	// Test each example separately
	for(std::vector<rna_example>::iterator example = rna_examples.begin(); example != rna_examples.end(); ++example)
	{
		Zuker zuker = Zuker(settings, (*example).sequence, thermodynamics);
		zuker.energy();
		zuker.traceback();
		
		std::string predicted_structure = "";
		zuker.dot_bracket.format((unsigned int)(*example).sequence.size() , predicted_structure);///@todo unsigned int -> size_t
		
		BOOST_REQUIRE_EQUAL((*example).dot_bracket_pattern.size() , predicted_structure.size());
		BOOST_CHECK_MESSAGE(db.match((*example).dot_bracket_pattern, predicted_structure) == false, "Predicted structure of '" << (*example).title << "' did match its true structure while it shouldn't:\n\t" << predicted_structure << " (predicted structure)\n\t" << (*example).dot_bracket_pattern << " (pattern of true structure)\n");
	}
}



/**
 * @brief Tests whether a specific sequence caused a critial error (seg-fault)
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_segfault_01)
{
	Sequence sequence = Sequence("CCCUUUGACCCAAAAGGGGCGAGGG");
	std::string true_structure = "(((...(((((......))))))))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	std::vector<rna_example> rna_examples;
	ReadSegments r = ReadSegments(settings.segment_filename);
	r.parse(thermodynamics.segments, thermodynamics.segmentloops, rna_examples);
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure); ///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_01a)
{
	Sequence sequence = Sequence("GGGAAAACCC");
	std::string true_structure = "(((....)))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -2.0999999f);
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_01b)
{
	Sequence sequence = Sequence("GGGAAAAACCC");
	std::string true_structure = "(((.....)))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -2.0999999f);
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_01c)
{
	Sequence sequence = Sequence("GGGAAAAAACCC");
	std::string true_structure = "(((......)))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -2.29999971f);//@mfold says -2.3
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_02)
{
	Sequence sequence = Sequence("GAAAAAAAAAAAAAAAAAAAAAAU");
	std::string true_structure = "........................";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , 0.0f);//@mfold says 0.0
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_03)
{
	Sequence sequence = Sequence("AGAAAAAAAAAAAAAAAAAAAAAAUA");
	std::string true_structure = "..........................";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , 0.0f);//@mfold says 0.0
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_04)
{
	Sequence sequence = Sequence("GGGGAAACCCCGGCGAAAAACGCC");// doesn't work while aGGGGAAACCCCAAAGGCGAAAAACGCCa does
	std::string true_structure = "((((...))))((((.....))))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -9.90000057f);
}




/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_05)
{
	Sequence sequence = Sequence("GGGGAAACCCCAGGCGAAAAACGCC");// doesn't work while aGGGGAAACCCCAAAGGCGAAAAACGCCa does
	std::string true_structure = "((((...)))).((((.....))))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -9.90000057f);
}




/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_06)
{
	Sequence sequence = Sequence("GGGGAAACCCCAAGGCGAAAAACGCC");// doesn't work while aGGGGAAACCCCAAAGGCGAAAAACGCCa does
	std::string true_structure = "((((...))))..((((.....))))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -9.90000057f);
}


/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_07)
{
	Sequence sequence = Sequence("GGGGAAACCCCAAAGGCGAAAAACGCC");// doesn't work while aGGGGAAACCCCAAAGGCGAAAAACGCCa does
	std::string true_structure = "((((...))))...((((.....))))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -9.90000057f);
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test_zuker_08)
{
	// Theoretically tRNA's can't be folded:
	// GGCCGGAAACCGGCCGCGCAAAAGCGCCCCGGGAAACCCGGG
	// But in practice it doesn't happen:
	Sequence sequence = Sequence("GCCGGCAAAGGCCGGAAACCGGCCAAGCGCAAAAGCGCAACCCGGGAAACCCGGGAAAGCCGGC");
	std::string true_structure = "((((((...((((((...))))))..((((....))))..((((((...))))))...))))))";
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -43.000004f);//-41.30 according to vienna
}



/**
 * @brief Tests prediction of certain structure tested with mfold/vienna
 *
 * @test
 * 
 * @todo figure out why the interior loop is chosen instead of the bulge - previous analysis indicated that the interior loop might get a too large energy score
 */
BOOST_AUTO_TEST_CASE(Test_zuker_09)
{
	Sequence sequence = Sequence("CCCaaaGGGGGAAACCCCCGGG");
	std::string true_structure = "(((...(((((...))))))))";
	//                        ->  .((...(((((...))))).))    << seems to be optimal - should be picked up as wij bifurcation because it is unpaired!
	
	// Load variables etc.
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	float energy = zuker.energy();
	
	// i = 2
	// i' = 6
	// j' = 18
	// j = 20
	
	zuker.traceback();
	std::string predicted_structure;
	zuker.dot_bracket.format((unsigned int) sequence.size() , predicted_structure);///@todo unsigned int -> size_t
	
	BOOST_CHECK_MESSAGE(predicted_structure.compare(true_structure) == 0, "Predicted structure '" << predicted_structure << "' and true structure '" << true_structure << "' are different");
	BOOST_CHECK_EQUAL(energy , -12.4000006f);
}



///@todo Test function with different minimum hairpin size -- should fail at the moment, because somewhere in the traceback or so the number 3 is hard coded

BOOST_AUTO_TEST_SUITE_END()
