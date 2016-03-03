# @file test/CMakeLists.txt
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
# Copyright (C) 2012-2016 Youri Hoogstrate
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

add_definitions(-std=c++0x)
add_definitions(-DBOOST_TEST_DYN_LINK)


set(BUILD_DIR      "../build")
set(BUILD_TEST_DIR "${BUILD_DIR}/test")


add_executable(direction_test          Direction_test.cpp)
add_executable(dotbracket_test         DotBracket_test.cpp ../src/DotBracket.cpp)
add_executable(nucleotide_test         Nucleotide_test.cpp)
add_executable(position_test           Position_test.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
add_executable(sequence_test           Sequence_test.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
add_executable(pairingplus_test        PairingPlus_test.cpp ../src/Sequence.cpp ../src/PairingPlus.cpp ../src/Pairing.cpp ../src/SubSequence.cpp)
add_executable(pairing_test            Pairing_test.cpp ../src/Pairing.cpp)
#add_executable(segmenttraceback_test   SegmentTraceback_test.cpp ../src/SegmentTraceback.cpp)
#add_executable(segment_test            Segment_test.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
#add_executable(segmentloop_test        SegmentLoop_test.cpp ../src/SegmentTraceback.cpp ../src/SegmentLoop.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
add_executable(segmenttreeelement_test SegmentTreeElement_test.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/Sequence.cpp  ../src/SubSequence.cpp)
add_executable(segmenttree_test        SegmentTree_test.cpp ../src/SegmentTree.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
add_executable(segmentlooptree_test    SegmentLoopTree_test.cpp ../src/SegmentTraceback.cpp ../src/SegmentLoop.cpp ../src/SegmentLoopTree.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
add_executable(readsegments_test       ReadSegments_test.cpp ../src/ReadSegments.cpp ../src/SegmentTree.cpp ../src/SegmentLoopTree.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/SegmentLoop.cpp ../src/Sequence.cpp ../src/PairingPlus.cpp ../src/Pairing.cpp ../src/SubSequence.cpp)
add_executable(settings_test           Settings_test.cpp ../src/Settings.cpp ../src/Utils/file_exists.cpp ../src/SegmentTree.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/Sequence.cpp ../src/SubSequence.cpp)
add_executable(scoringmatrix_test      ScoringMatrix_test.cpp ../src/ScoringMatrix.cpp)
add_executable(gibbsfreeenergy_test    GibbsFreeEnergy_test.cpp ../src/GibbsFreeEnergy.cpp ../src/ReadData.cpp ../src/SegmentTree.cpp ../src/SegmentLoopTree.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/SegmentLoop.cpp ../src/Sequence.cpp ../src/Pairing.cpp ../src/PairingPlus.cpp ../src/SubSequence.cpp )
add_executable(zuker_test              Zuker_test.cpp ../src/Zuker.cpp ../src/GibbsFreeEnergy.cpp ../src/ReadData.cpp ../src/ReadSegments.cpp ../src/SegmentTree.cpp ../src/SegmentLoopTree.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/SegmentLoop.cpp ../src/Sequence.cpp ../src/Pairing.cpp ../src/PairingPlus.cpp ../src/SubSequence.cpp ../src/Settings.cpp ../src/Utils/file_exists.cpp ../src/ScoringMatrix.cpp ../src/DotBracket.cpp)
add_executable(zuker_traceback_test    Zuker_traceback_test.cpp ../src/Zuker.cpp ../src/GibbsFreeEnergy.cpp ../src/ReadData.cpp ../src/ReadSegments.cpp ../src/SegmentTree.cpp ../src/SegmentLoopTree.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/SegmentLoop.cpp  ../src/Sequence.cpp ../src/Pairing.cpp ../src/PairingPlus.cpp ../src/SubSequence.cpp ../src/Settings.cpp ../src/Utils/file_exists.cpp ../src/ScoringMatrix.cpp ../src/DotBracket.cpp)
add_executable(main_test               main_test.cpp ../src/Zuker.cpp ../src/GibbsFreeEnergy.cpp ../src/ReadData.cpp ../src/ReadSegments.cpp ../src/SegmentTree.cpp ../src/SegmentLoopTree.cpp ../src/SegmentTreeElement.cpp ../src/SegmentTraceback.cpp ../src/Segment.cpp ../src/SegmentLoop.cpp ../src/Sequence.cpp ../src/Pairing.cpp ../src/PairingPlus.cpp ../src/SubSequence.cpp ../src/Settings.cpp ../src/Utils/file_exists.cpp ../src/ScoringMatrix.cpp ../src/DotBracket.cpp)


set_target_properties(direction_test
                      dotbracket_test
                      nucleotide_test
                      position_test
                      sequence_test
                      pairingplus_test
                      pairing_test
#                      segmenttraceback_test
#                      segment_test
#                      segmentloop_test
                      segmenttreeelement_test
                      segmenttree_test
                      segmentlooptree_test
                      readsegments_test
                      settings_test
                      scoringmatrix_test
                      gibbsfreeenergy_test
                      zuker_test
                      zuker_traceback_test
                      main_test
                                      PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${BUILD_TEST_DIR}")
