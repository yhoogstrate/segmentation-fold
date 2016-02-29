/**
 * @file src/SegmentTraceback.cpp
 *
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



#include "main.hpp"

#include <array>

#include "Pair.hpp"
#include "SegmentTraceback.hpp"



/**
 * @brief Constructor of the Segment class.
 *
 * @section DESCRIPTION
 * <PRE>
 * 5') ............(i)...(i')..
 *       |||  |||   | ::   |   .
 * 3') ............(j)...(j')..
 * </PRE>
 *
 * tb:  --------->  [i,j], [i+1, j-1], [i+2, j-2] -----> [i', j']...
 *
 *
 * Double check the following:
 *
 *     Alignment:       Bonds from i,j:    Traceback from i,j:
 * <PRE>
 * 5') (i) UGUGAUA      3,11               4,(-)3
 *            ::: A     4,10               5,(-)2
 * 3') (j)    AGUA      5,9                6,(-)1
 * </PRE>
 *
 * <PRE>
 * 5') ............(i)...
 *       |||  |||   | :: .
 * 3') ............(j)...
 *</PRE>
 * tb:  --------->  [i,j], [i+1, j-1], [i+2, j-2]
 *
 * @param arg_bonds The bonds are given from outside to inside (FiFo)
 *
 */
SegmentTraceback::SegmentTraceback(std::vector <Pair> arg_bonds):
	bonds(arg_bonds)
{
	this->reset();
}



/**
 * @brief Resets the traceback to the first bond of the segment
 *
 */
void SegmentTraceback::reset(void)
{
	this->it = this->bonds.begin();
}



/**
 * @brief Returns the number of bonds present in the traceback
 *
 * */
size_t SegmentTraceback::size(void)
{
	return this->bonds.size();
}



/**
 * @brief FiFo popping of the Pairs, + updating iterators i and j and reset the tracebac iterator if the end has reached
 *
 * @param i Reference to variable i that should be traced back
 * @param j Reference to variable j that should be traced back
 *
 * @return Whether the return was VALID - if false is returned, the traceback is being reset to it's origin and will return true again.
 *
 */
bool SegmentTraceback::traceback(unsigned int &i, unsigned int &j)
{
	if(this->it == this->bonds.end())
	{
		this->reset();
		
		return false;
	}
	else
	{
		i += (*this->it).first;
		j -= (*this->it).second;
		
		this->it++;
		
		return true;
	}
}
