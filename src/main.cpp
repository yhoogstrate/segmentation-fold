/**
 * @file src/main.cpp
 *
 * @section DESCRIPTION
 * This file contains the main() function.
 * The program implements the Zuker's minimum free energy model with
 * Segment/K-turn functionality for RNA secondary structure prediction.
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2016 Youri Hoogstrate
 *
 * This file is part of segmentation-fold
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

#include "Settings.hpp"
#include "DotBracket.hpp"
#include "ReadData.hpp"
#include "ScoringMatrix.hpp"
#include "Zuker.hpp"



/**
 * @brief Predicts the structure of an RNA sequence in major 3 steps: 1) obtain parameters from IO, 2) Run the Zuker algorithm, 3) print the 2D structure.
 *
 * @param argc Number of commandline arguments.
 * @param argv Array of strings with the actual commandline arguments.
 *
 * @return 0 For success; 1 for abnormal termination
 */
int main(int argc, char *argv[])
{
	// Load data
	Sequence sequence = Sequence();
	
	// Set variables
	Settings settings = Settings(argc, argv, sequence);
	if(!settings.run_print_version and !settings.run_print_usage)
	{
		ReadData thermodynamics = ReadData();
		ReadSegments readsegments = ReadSegments(settings.segment_filename);
		
		if(settings.segment_prediction_functionality)
		{
			readsegments.parse(thermodynamics.segments, thermodynamics.segmentloops);
		}
		
		// Run algorithm
		Zuker zuker = Zuker(settings, sequence, thermodynamics);		// Zuker algorithm
		zuker.energy();													// - filling phase
		zuker.traceback();												// - traceback
		zuker.print_2D_structure();
	}
	
	return 0;
}
