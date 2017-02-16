/**
 * @file include/SegmentTraceback.hpp
 *
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
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
 */


#ifndef SEGMENTTRACEBACK_HPP
#define SEGMENTTRACEBACK_HPP


#include <vector>


/**
 * @brief Container for traceback coordinates for segment objects, to avoid pathcrossing
 *
 * @section DESCRIPTION
 * A container of coordinates that allow traceback with non-canonical
 * pairs, derived from either a Segment or SegmentLoop object.
 *
 * The first implementation calculated back from [i',j'] but since a
 * loop only has [i,j] and not [i',j'], this traceback should be
 * calculating back relative to [i,j].
 *
 * ---- old style, relative to [i',j'] ----
 * Alignment:    Bonds from i j :    Traceback from i',j':
 * 5') UGUGAU    3-2                 -1,-1
 *        |||    4-1                 -2,-2
 * 3')    AGU    5-0                 -3,-3
 *
 * ---- new style, relative to [i,j] ----
 * Alignment:    Bonds from i j*:    Traceback from i,j:
 * 5') UGUGAU    3-8                 (i+)4,(j-)1
 *        |||    4-7                 (i+)5,(j-)2
 * 3')    AGU    5-6                 (i+)6,(j-)3
 *
 * *pretend that the sequence is 5') UGUGAU UGA (3'
 *
 *
 * 5') ............(i)...(i')..
 *       |||  |||   | ::   |   .
 * 3') ............(j)...(j')..
 *
 * tb:  --------->  [i,j], [i+1, j-1], [i+2, j-2] -----> [i', j']...
 *
 *
 *
 * Double check the following:
 *
 *     Alignment:       Bonds from i,j:    Traceback from i,j:
 * 5') (i) UGUGAUA      3,11               4,(-)3
 *            ::: A     4,10               5,(-)2
 * 3') (j)    AGUA      5,9                6,(-)1
 *
 * 5') ............(i)...
 *       |||  |||   | :: .
 * 3') ............(j)...
 *
 * tb:  --------->  [i,j], [i+1, j-1], [i+2, j-2]
 *
 *
 */
class SegmentTraceback
{
		friend class Test_SegmentTraceback;
		
	private:
		std::vector<Pair>::iterator it;
		std::vector<Pair> bonds;
		
	public:
		SegmentTraceback(std::vector<Pair> arg_bonds);
		
		bool traceback(unsigned int &i, unsigned int &j);
		void reset(void);
		
		size_t size(void);
};



///@brief Friend class of SegmentTraceback that allows testing its private members
class Test_SegmentTraceback: public SegmentTraceback
{
	public:
		using SegmentTraceback::SegmentTraceback;
		using SegmentTraceback::it;
		using SegmentTraceback::bonds;
};



#endif	// SEGMENTTRACEBACK_HPP
