/**
 * @file src/SegmentLoopTree.cpp
 *
 * @date 2015-12-05
 *
 * @author Youri Hoogstrate
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



#include "main.hpp"
#include "Utils/utils.hpp"

#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Position.hpp"
#include "SubSequence.hpp"
#include "Sequence.hpp"

#include "SegmentLoop.hpp"
#include "SegmentLoopTree.hpp"



/**
 * @brief Initializes an empty tree
 *
 * @date 2015-12-05
 */
SegmentLoopTree::SegmentLoopTree()
{
	this->root = nullptr;
}



/**
 * @brief Destructs the tree recursively
 *
 * @date 2015-12-05
 */
SegmentLoopTree::~SegmentLoopTree()
{
	delete this->root;
}



/**
 * @brief Searches a pair by using sets of iterators
 *
 * @date 2015-12-05
 *//*
Segment *SegmentLoopTree::search(SubSequence &arg_segment5p, SubSequence &arg_segment3p)
{
	return (this->empty()) ? nullptr : this->root->search_segment(arg_segment5p, arg_segment3p);
}*/



/**
 * @brief Adds a Segment (as member of a SegmentLoopTreeElement) to the tree.
 *
 * @date 2015-12-05
 */
void SegmentLoopTree::insert(Segment &arg_segment)
{
	if(this->empty())
	{
		SegmentLoopTreeElement *element;
		element = new SegmentLoopTreeElement(arg_segment, nullptr, nullptr);
		this->root = element;
	}
	else
	{
		this->insert(arg_segment,this->root);
	}
}



/**
 * 
 * @date 2015-12-05
 */
void SegmentLoopTree::insert(SegmentLoop &arg_segmentloop, SegmentLoopTreeElement *arg_element)
{
	switch ( (arg_element->sequence.compare(arg_segmentloop.sequence))
	{
		case IS_SMALLER:
			if(arg_element->left!=nullptr)
			{
				insert(arg_segmentloop, arg_element->left);
			}
			else
			{
				arg_element->left=new SegmentLoopTreeElement(arg_segmentloop, nullptr, nullptr);
			}
			break;
		
		case IS_LARGER:
			if(arg_element->right!=nullptr)
			{
				insert(arg_segmentloop, arg_element->right);
			else
			{
				arg_element->right=new SegmentLoopTreeElement(arg_segmentloop, nullptr, nullptr);
			}
		break;
		
		//default:  inserting same element twice
		// if debug: throw exception
	}
}


/**
 * @brief Returns whether the tree is empty or not
 *
 * @date 2015-05-02
 */
bool SegmentLoopTree::empty(void)
{
	return (this->root == nullptr);
}



/**
 * @brief Counts recursively the elements in the tree
 *
 * @date 2015-12-05
 *
 * @todo If this function becomes inline, it can not be found anymore... This is probably because its kinda recursive since it uses another inline function? Inlines seem to have to be defined in a header file
 * @link http://www.parashift.com/c++-faq-lite/inline-member-fns.html
 */
size_t SegmentLoopTree::size(void)
{
	return this->empty() ? 0 : this->root->size();
}
