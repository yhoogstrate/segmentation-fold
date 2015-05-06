/**
 * @file include/SegmentTreeElement.hpp
 *
 * @date 2015-05-02
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



#ifndef SEGMENTTREEELEMENT_HPP
#define	SEGMENTTREEELEMENT_HPP

#include "SubSequence.hpp"

/**
 * @brief Tree element for fast searching, and contains a reference to the segment object.
 *
 * @date 2015-04-29
 */
class SegmentTreeElement
{
	private:
		Segment &segment;
		
		SegmentTreeElement *link_5p_size_smaller;
		SegmentTreeElement *link_5p_size_larger;
		
		SegmentTreeElement *link_5p_segment_smaller;
		SegmentTreeElement *link_5p_segment_larger;
		
		
		SegmentTreeElement *link_3p_size_smaller;
		SegmentTreeElement *link_3p_size_larger;
		
		SegmentTreeElement *link_3p_segment_smaller;
		SegmentTreeElement *link_3p_segment_larger;
		
		void add_segment(Segment &arg_segment, char &arg_search_type);
		Segment *search_segment(SubSequence &arg_sequence_5p, SubSequence &arg_sequence_3p, char &arg_search_type);
		
	public:
		SegmentTreeElement(Segment &arg_segment);
		~SegmentTreeElement(void);
		
		void add_segment(Segment &arg_segment);
		Segment *search_segment(SubSequence &arg_sequence_5p, SubSequence &arg_sequence_3p);
		
		size_t size(void);
};


#endif	// SEGMENTTREEELEMENT_HPP
