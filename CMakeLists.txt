# @file CMakeLists.txt
#
# @date 15-apr-2015
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
file(MAKE_DIRECTORY "${BUILD_DIR}")
file(MAKE_DIRECTORY "${BUILD_DIR}/test")

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
# ------------------------------ Styling -------------------------------


# ----------------------------------------------------------------------
# --------------------------- Documentation ----------------------------


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

add_test(direction_test          build/test/direction_test)
# 
add_test(nucleotide_test         build/test/nucleotide_test)
add_test(position_test           build/test/position_test)
add_test(sequence_test           build/test/sequence_test)
add_test(segment_test            build/test/segment_test)

# ----------------------------------------------------------------------
# --------------------------- Installation -----------------------------

# ----------------------------------------------------------------------
