/**
 * @file include/ReadData.hpp
 *
 * @date 2015-12-07
 *
 * @author Youri Hoogstrate
 * @author Lisa Yu
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2015 Youri Hoogstrate
 *
 * This file is part of segmentation-fold and originally taken from
 * yh-kt-fold.
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



#ifndef READDATA_HPP
#define	READDATA_HPP



#include "main.hpp"

#include "../include/Utils/utils.hpp"

#include <vector>
#include <algorithm>
#include <map>


/**
 * @brief Reads the thermodynamic Gibbs energy parameters.
 *
 * @date 2015-12-07
 *
 * @section DESCRIPTION
 * Read thermodynamics data (Turners dataset for now, Andronescu maybe
 * in the future?)
 *
 * The following variables are enlisted:
 * - bulge			Energy-values for bulge loops
 * - dangle
 * - miscloop
 * - hairpin		Energy-values for hairpin loops
 * - inter			Energy-values for interior/internal loops
 * - poppen
 *
 * - int11			Energy-values for Symetrical interior loops	(2 stands probably for number of pairing surrounding the interior loop?)
 * - asint3			Energy-values for A-symetrical interior loops (odd number of mismatches??)
 * - int22			Energy-values for Symetrical interior loops	(4 stands probably for number of pairing surrounding the interior loop?)
 * - asint5			Energy-values for A-symetrical interior loops (odd number of mismatches??)
 * - sint6			Energy-values for Symetrical interior loops	(6 stands probably for number of pairing surrounding the interior loop?)
 *
 * - stack
 *
 * - tstkh
 * - tstki
 *
 * - tloop
 * - triloop
 *
 * - prelog
 *
 * - segments			List of all K-turns found in the data-file
 */
class ReadData
{
	public:
		float int11[6][6][4][4],										// formerly known as sint4
			  int21[6][6][4][4][4],										// formerly known as asint5
			  int22[6][6][4][4][4][4],									// formerly known as sint6
			  
			  stack[6][6],												// Stack that closes an other stack
			  tstackh[6][4][4],											// Stack that closes a hairpin
			  tstki[6][4][4],											// Stack that closes an interior loop
			  
			  miscloop[20];
			  
		std::map<Sequence, float> tloop_map;
		std::map<Sequence, float> triloop_map;
		
		std::vector<float> loop_hairpin;
		std::vector<float> loop_hairpin_C_penalty;
		std::vector<float> loop_bulge;
		std::vector<float> loop_interior;
		
		std::vector<float> poppen_p;
		
		SegmentTree segments;
		SegmentLoopTree segmentloops;
		
		unsigned int minimal_hairpin_length;							///@todo this is an algorithmic setting - move it to Settings
		
		
		ReadData();
		
		char **get_datadirs(void);
		std::string *get_segments_file(void);
		
		void load_int11();
		void load_int21();
		void load_int22();
		
		void load_tstackh();
		void load_tstacki();
		
		void load_loop_hairpin();
		void load_loop_interior();
		void load_loop_bulge();
		
		void load_poppen();
		void load_stack();
		void load_triloop();
		void load_tloop();
		void load_miscloop();
};

#endif	// READDATA_HPP
