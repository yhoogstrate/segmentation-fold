# @file CMakeLists.txt
#
# @date 17-apr-2015
#
# @author Youri Hoogstrate
#
# @section DESCRIPTION
# This CMakefile contains the description of compilance of the tool
# segmentation-fold as well as building several unit- # and functional
# tests. By running 'cmake .' followed by the 'make check' command, the
# software can simply be checked for (new) introduced mistakes.
#
# @section LICENSE
# segmentation-fold can predict RNA 2D structures including K-turns.
# Copyright (C) 2012-2015 Youri Hoogstrate
#
# This file is part of segmentation-fold.
#
# segmentation-fold is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# segmentation-fold is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

# ----------------------------- Packaging ------------------------------

project(segmentation-fold)

set(PROJECT_VERSION "1.0.0")
set(PACKAGE_URL "https://github.com/yhoogstrate/segmentation-fold")
set(PACKAGE_BUGREPORT "${PACKAGE_URL}/issues")

set(MOTIFS_FILE "motifs.xml")

set(BUILD_DIR "build")
set(BUILD_TEST_DIR "${BUILD_DIR}/test")
file(MAKE_DIRECTORY "${BUILD_DIR}")
file(MAKE_DIRECTORY "${BUILD_TEST_DIR}")

# include the git commit SHA1
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
if(NOT "${GIT_SHA1}" STREQUAL "")
	set(GIT_SHA1_STRING "-${GIT_SHA1}")
endif()

# ----------------------------------------------------------------------
# ------------------------------ Pre-comp ------------------------------

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE debug)
	message(WARNING "No CMAKE_BUILD_TYPE was defined. Build type 'debug' was chosen, which will perform slower than using 'cmake -DCMAKE_BUILD_TYPE=release'.")
endif(NOT CMAKE_BUILD_TYPE)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE STREQUAL "debug")
	set(DEBUG "true")
else()
	set(DEBUG "false")
endif()

configure_file("include/config.hpp.in" "include/config.hpp")

# ----------------------------------------------------------------------
# ------------------------------ README --------------------------------
# run 'make readme' to create/update the README.md file

add_custom_command(
				OUTPUT  make_readme
				
				COMMAND cat doc/md/INTRODUCTION.md doc/md/INSTALL.md > README.tmp
				COMMAND mv README.tmp README.md
				
				COMMAND cat README.md doc/md/USAGE.md > README.tmp
				COMMAND mv README.tmp README.md
				
				COMMAND cat README.md AUTHORS > README.tmp
				COMMAND mv README.tmp README.md
				
				COMMAND cat README.md CONTRIBUTING.md > README.tmp
				COMMAND mv README.tmp README.md
			)

add_custom_target(readme DEPENDS make_readme)

# ----------------------------------------------------------------------
# ------------------------------ Styling -------------------------------
# run 'make tidy' to style the syntax using astyle

#find_program(ASTYLE astyle)

add_custom_command(
				OUTPUT  make_tidy
				COMMAND astyle --options=share/.astylerc include/*.hpp
				COMMAND astyle --options=share/.astylerc src/*.cpp src/Utils/*.cpp
				COMMAND astyle --options=share/.astylerc test/*.cpp
				COMMAND date > build/astyle.timestamp.txt
			)

add_custom_target(tidy DEPENDS make_tidy )

# ----------------------------------------------------------------------
# --------------------------- Documentation ----------------------------

#find_program(ASTYLE doxygen)

configure_file("doc/Doxyfile.in" "doc/Doxyfile")

add_custom_command(
				OUTPUT  make_doc
				COMMAND doxygen doc/Doxyfile
			)

add_custom_target(doc DEPENDS make_doc)

# ----------------------------------------------------------------------
# ---------------------------- Compilation -----------------------------

add_subdirectory(src)
include_directories(include)
add_definitions(-std=c++11)

# ----------------------------------------------------------------------
# ------------------------------ Testing -------------------------------

add_subdirectory(test)

enable_testing()

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})					# 'make check' as alias for 'make test'

add_test(direction_test          "${BUILD_TEST_DIR}/direction_test")
add_test(nucleotide_test         "${BUILD_TEST_DIR}/nucleotide_test")
add_test(position_test           "${BUILD_TEST_DIR}/position_test")
add_test(sequence_test           "${BUILD_TEST_DIR}/sequence_test")
add_test(segment_test            "${BUILD_TEST_DIR}/segment_test")

# ----------------------------------------------------------------------
# --------------------------- Installation -----------------------------

# ----------------------------------------------------------------------
