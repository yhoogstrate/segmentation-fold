/**
 * @file test/Settings_test.cpp
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



#define BOOST_TEST_MODULE Settings


#include <stdio.h>
#include <stdlib.h>



#include "Utils/utils.hpp"

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

#include "ReadData.hpp"
#include "Settings.hpp"

#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether a sequence can be obtained from argv correctly (-s)
 *
 * @test
 *
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "ACTGactgACUGacug", nullptr};
	signed int argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
	
	Sequence sequence;
	
	Settings settings = Settings(argc, argv, sequence);
	
	BOOST_REQUIRE_EQUAL(sequence.size(), 16);
	
	size_t i = 0;
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK(sequence == Sequence("ACTGactgACUGacug"));
}

/**
 * @brief Tests whether a FASTA file can be parsed correctly (-f)
 *
 * @test
 *
 *
 * @todo Implement the possibility to run multiple entries from a FASTA file
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	std::string filename = "tmp.settings_test_test2";
	
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile <<   "\n\n  \n>some_sequence\nACTG\nactg\nACUG\nacug\n\n>next seq\nacgtgactgac\n";
	myfile.close();
	
	
	char *argv[] = {(char *) PACKAGE_NAME, (char *) "-f", (char *) filename.c_str(), nullptr};
	signed int argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
	
	Sequence sequence;
	
	Settings settings = Settings(argc, argv, sequence);
	
	BOOST_REQUIRE_EQUAL(sequence.size(), 16);
	
	size_t i = 0;
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::A);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::C);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::U);
	BOOST_CHECK_EQUAL(sequence[i++], Nucleotide::G);
	
	BOOST_CHECK(sequence == Sequence("ACTGactgACUGacug"));
	
	unlink(filename.c_str());
}

/**
 * @brief Tests whether segment functionality can be en/disabled (-p)
 *
 * @test
 *
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	Sequence sequence;
	signed int argc;
	std::string is;
	
	{
		// Check default value
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.segment_prediction_functionality, true);
	}
	
	{
		// Check enabling
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", (char *) "-p", (char *) "1", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.segment_prediction_functionality, true);
	}
	
	{
		// Check disabling
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", (char *) "-p", (char *) "0", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.segment_prediction_functionality, false);
	}
}



/**
 * @brief Tests whether the minimal hairpin-size can be set (-H)
 *
 * @test
 *
 */
BOOST_AUTO_TEST_CASE(Test4)
{
	Sequence sequence;
	int argc;
	std::string is;
	
	
	// Check default value
	char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", nullptr};
	argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
	
	Settings settings = Settings(argc, argv, sequence);
	
	BOOST_CHECK_EQUAL(settings.minimal_hairpin_length, 3);
	
	
	// Check argumented values
	for(unsigned int i = 0; i < 100; i ++)
	{
		is = std::to_string(i);
		char *ics = (char *) is.c_str();
		
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", (char *) "-H", (char *) ics, nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_MESSAGE(settings.minimal_hairpin_length == i, "Failed to obtain min_hairpin_size of " << i);
	}
}

/**
 * @brief Tests whether correct segment xml file is chosen (-x)
 *
 * @test
 *
 */
BOOST_AUTO_TEST_CASE(Test5)
{
	Sequence sequence;
	int argc;
	
	{
		// Check example file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", (char *) "-x", (char *) "share/segmentation-fold/" SEGMENTS_FILE, nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.segment_filename, "share/segmentation-fold/" SEGMENTS_FILE);
	}
	
	{
		// Check non existing file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", (char *) "-x", (char *) "/dev/null/neverexist", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		BOOST_CHECK_THROW(Settings settings = Settings(argc, argv, sequence), std::invalid_argument);
	}
}

/**
 * @brief Tests whether the version command is being picked up (-V)
 *
 * @test
 *
 */
BOOST_AUTO_TEST_CASE(Test6)
{
	Sequence sequence;
	int argc;
	
	{
		// Check example file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, true);
	}
	
	{
		// Check non existing file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", (char *) "-V", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, false);
	}
}

/**
 * @brief Tests whether the version command is being picked up (--version)
 *
 * @test
 *
 */
BOOST_AUTO_TEST_CASE(Test7)
{
	Sequence sequence;
	int argc;
	
	{
		// Check example file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, true);
	}
	
	{
		// Check non existing file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "--version", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, false);
	}
}

/**
 * @brief Tests whether the version command is being picked up (--help)
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test8)
{
	Sequence sequence;
	int argc;
	
	{
		// Check example file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, true);
	}
	
	{
		// Check non existing file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "--help", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, false);
	}
	
	{
		// Check non existing file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-h", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, false);
	}
}


/**
 * @brief Tests whether the number of threads is being picked up (-t)
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test9)
{
	Sequence sequence;
	int argc;
	std::string is;
	
	
	// Check default value
	char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", nullptr};
	argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
	
	Settings settings = Settings(argc, argv, sequence);
	
	BOOST_CHECK_EQUAL(settings.num_threads, 0);
	
	
	// Check argumented values
	for(unsigned int i = 0; i < 100; i ++)
	{
		is = std::to_string(i);
		char *ics = (char *) is.c_str();
		
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", (char *) "-t", (char *) ics, nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_MESSAGE(settings.num_threads  == i, "Failed to obtain num_threads of " << i << " ( " << settings.num_threads << " was found instead)");
	}
}


/**
 * @brief Tests whether the default xml command is being picked up (-X, --default-xml)
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test10)
{
	Sequence sequence;
	int argc;
	
	{
		// Check example file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-s", (char *) "a", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, true);
	}
	
	{
		// Check non existing file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "--default-xml", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, false);
	}
	
	{
		// Check non existing file
		char *argv[] = {(char *) PACKAGE_NAME, (char *) "-X", nullptr};
		argc = (signed int) sizeof(argv) / (signed int) sizeof(char *) - 1;
		
		Settings settings = Settings(argc, argv, sequence);
		
		BOOST_CHECK_EQUAL(settings.proceed_with_folding, false);
	}
}

BOOST_AUTO_TEST_SUITE_END()
