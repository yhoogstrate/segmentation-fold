/**
 * @file include/SegmentLoopTree.hpp
 *
 * @date 2015-12-05
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


#ifndef SEGMENTTREE_HPP
#define	SEGMENTTREE_HPP


#include "SubSequence.hpp"



struct SegmentLoopTreeElement
{
	SegmentLoop &segmentloop;
	
	SegmentLoopTreeElement *left;// = nullptr;
	SegmentLoopTreeElement *right;// = nullptr;
};



/**
 * @brief Stores and provides access to (K-loop) segmentloop objects and provides extended trace-back routes for such segments.
 *
 * @section DESCRIPTION
 * Stores all (K-turn) segments in memory.
 *
 * @date 2015-12-05
 */
class SegmentLoopTree
{
	private:
		SegmentLoopTreeElement *root;
		
	public:
		SegmentLoopTree();
		~SegmentLoopTree();
		
		void insert(SegmentLoop &arg_segment);
		void insert(SegmentLoop &arg_segment, SegmentLoopTreeElement *arg_element);
		
		SegmentLoop *search(SubSequence &arg_segment);
		
		bool empty(void);
		size_t size(void);
};


#endif	// SEGMENTTREE_HPP
