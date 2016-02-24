/**
 * @file test/DotBracket_test.cpp
 *
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



#define BOOST_TEST_MODULE DotBracket



#include "DotBracket.hpp"
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_SUITE(Testing)

/**
 * @brief Tests whether the match function matches an exactly identical string
 *
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	DotBracket db = DotBracket();
	
	std::string dot_bracket_pattern = "(((...)))";
	std::string dot_bracket_subject = "(((...)))";
	
	BOOST_CHECK_MESSAGE(db.match(dot_bracket_pattern, dot_bracket_subject), "DotBracket pattern '" << dot_bracket_pattern << "' doesn't match subject '" << dot_bracket_subject << "'");
}

/**
 * @brief Tests whether the match function matches two different strings as different
 *
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	DotBracket db = DotBracket();
	
	std::string dot_bracket_pattern = "(((...())))";
	std::string dot_bracket_subject = "(((()...)))";
	
	BOOST_CHECK_MESSAGE(db.match(dot_bracket_pattern, dot_bracket_subject) == false, "DotBracket pattern '" << dot_bracket_pattern << "' does match subject '" << dot_bracket_subject << "'");
}


/**
 * @brief Tests whether the match function matches a structure that matches a pattern
 *
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	DotBracket db = DotBracket();
	
	std::string dot_bracket_pattern = ".(((\?\?\?\?\?\?\?\?\?\?)))..()";///@note Question marks are neccesairy because they will be interpreted as some binary flags otherwise
	std::string dot_bracket_subject = ".(((()..(..()))))..()";
	
	BOOST_CHECK_MESSAGE(db.match(dot_bracket_pattern, dot_bracket_subject), "DotBracket pattern '" << dot_bracket_pattern << "' doesn't match subject '" << dot_bracket_subject << "'");
}

BOOST_AUTO_TEST_SUITE_END()
