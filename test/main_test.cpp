/**
 * @file test/main_test.cpp
 *
 * @date 2015-06-28
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

/**
 * @brief Tests whether an unstructured RNA is indeed unstructured
 *
 * @test
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
	zuker.dot_bracket.format(sequence.size() , predicted_structure);
	
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
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format(sequence.size() , predicted_structure);
	
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
 */
BOOST_AUTO_TEST_CASE(Test_bulge_loop)
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
 */
BOOST_AUTO_TEST_CASE(Test_interior_loop)
{
	// Initialize variables
	Sequence sequence = Sequence("GGGGAAAGGGGAAACCCCAAACCCC");
	std::string true_structure = "((((...((((...))))...))))";
	
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format(sequence.size() , predicted_structure);
	
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
 */
BOOST_AUTO_TEST_CASE(Test_bifurcation)
{
	// Initialize variables
	Sequence sequence = Sequence("GGGGgggaaaCCCgggAAAcccCCCC");
	std::string true_structure = "(((((((...)))(((...)))))))";
	
	Settings settings = Settings(0, nullptr, sequence);
	ReadData thermodynamics = ReadData();
	
	// Predict structure
	Zuker zuker = Zuker(settings, sequence, thermodynamics);
	zuker.energy();
	zuker.traceback();
	
	// Obtain and compare results
	std::string predicted_structure;
	zuker.dot_bracket.format(sequence.size() , predicted_structure);
	
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
	ReadSegments r = ReadSegments(settings.segment_filename, thermodynamics.segments, rna_examples);
	
	DotBracket db = DotBracket();
	
	// Test each example separately
	for(std::vector<rna_example>::iterator example = rna_examples.begin(); example != rna_examples.end(); ++example)
	{
		Zuker zuker = Zuker(settings, (*example).sequence, thermodynamics);
		zuker.energy();
		zuker.traceback();
		
		std::string predicted_structure = "";
		zuker.dot_bracket.format((*example).sequence.size() , predicted_structure);
		
		BOOST_REQUIRE_EQUAL((*example).dot_bracket_pattern.size() , predicted_structure.size());
		BOOST_CHECK_MESSAGE(db.match((*example).dot_bracket_pattern, predicted_structure), "Predicted structure of '" << (*example).title << "' doesn't match it's true structure:\n\t[" << predicted_structure << "] (predicted structure)\n\t[" << (*example).dot_bracket_pattern << "] (pattern of true structure)\n");
	}
}


///@todo Test function with segment-functionality disabled
///@todo Test function with different minimum hairpin size

BOOST_AUTO_TEST_SUITE_END()