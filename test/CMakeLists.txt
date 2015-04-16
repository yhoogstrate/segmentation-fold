# @file test/CMakeLists.txt
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

include_directories(../include)

find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)



add_executable(../build/test/direction_test Direction_test.cpp)
# 
add_executable(../build/test/nucleotide_test Nucleotide_test.cpp)
add_executable(../build/test/position_test Position_test.cpp ../src/Sequence.cpp)
add_executable(../build/test/sequence_test Sequence_test.cpp ../src/Sequence.cpp)
add_executable(../build/test/segment_test Segment_test.cpp ../src/Segment.cpp ../src/Sequence.cpp)
