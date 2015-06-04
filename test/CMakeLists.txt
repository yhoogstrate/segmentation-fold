# @file test/CMakeLists.txt
#
# @date 2015-06-04
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
add_definitions(-std=c++11)

find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)


set(BUILD_DIR      "../build")
set(BUILD_TEST_DIR "${BUILD_DIR}/test")


add_executable(direction_test          Direction_test.cpp)
add_executable(dotbracket_test         DotBracket_test.cpp ../src/DotBracket.cpp)
add_executable(nucleotide_test         Nucleotide_test.cpp)
add_executable(position_test           Position_test.cpp ../src/Sequence.cpp)
add_executable(sequence_test           Sequence_test.cpp ../src/Sequence.cpp)
add_executable(pairingplus_test        PairingPlus_test.cpp ../src/Sequence.cpp ../src/PairingPlus.cpp)
add_executable(segment_test            Segment_test.cpp ../src/Segment.cpp ../src/Sequence.cpp)
add_executable(segmenttreeelement_test SegmentTreeElement_test.cpp ../src/SegmentTreeElement.cpp ../src/Segment.cpp ../src/Sequence.cpp  ../src/SubSequence.cpp)
add_executable(segmenttree_test        SegmentTree_test.cpp ../src/SegmentTree.cpp ../src/SegmentTreeElement.cpp ../src/Segment.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
add_executable(readsegments_test       ReadSegments_test.cpp ../src/ReadSegments.cpp ../src/SegmentTree.cpp ../src/SegmentTreeElement.cpp ../src/Segment.cpp ../src/Sequence.cpp ../src/PairingPlus.cpp ../src/SubSequence.cpp)

set_target_properties(direction_test
                      dotbracket_test
                      nucleotide_test
                      position_test
                      sequence_test
                      pairingplus_test
                      segment_test
                      segmenttreeelement_test
                      segmenttree_test
                      readsegments_test
                                      PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${BUILD_TEST_DIR}")
