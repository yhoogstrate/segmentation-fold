/**
 * @file src/SegmentLoopTree.cpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2016 Youri Hoogstrate
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
 */
SegmentLoopTree::SegmentLoopTree()
{
	this->root = nullptr;
}



/**
 * @brief Destructs the tree recursively
 */
SegmentLoopTree::~SegmentLoopTree()
{
	this->clear();
}



/**
 * @brief Searches for a certain SegmentLoop given a SubSequence (2 pointers make it much quicker than a deep copy of a Sequence)
 */
SegmentLoop *SegmentLoopTree::search(SubSequence &arg_subsequence)
{
	return this->search(arg_subsequence, this->root);
}



/**
 * @brief Search per element
 */
SegmentLoop *SegmentLoopTree::search(SubSequence &arg_subsequence, SegmentLoopTreeElement *arg_element)
{
	if(arg_element != nullptr)
	{
		switch(arg_element->segmentloop.sequence.compare(arg_subsequence))
		{
			case IS_SMALLER:// tree element is smaller, thus the subsequence is right (larger) of the current element
				return this->search(arg_subsequence, arg_element->right);
				break;
			case IS_LARGER:// the tree element is larger, thus the subsequence is left (smaller) of the current element
				return this->search(arg_subsequence, arg_element->left);
				break;
			default://case IS_EQUAL - return current SegmentLoop:
				return &arg_element->segmentloop;
				break;
		}
	}
	else
	{
		return nullptr;
	}
}



/**
 * @brief Adds a Segment (as member of a SegmentLoopTreeElement) to the tree.
 */
void SegmentLoopTree::insert(SegmentLoop &arg_segmentloop)
{
	if(this->empty())
	{
		SegmentLoopTreeElement *element;
		element = new SegmentLoopTreeElement(arg_segmentloop);
		this->root = element;
	}
	else
	{
		this->insert(arg_segmentloop, this->root);
	}
}



/**
 * @brief Insert element
 */
void SegmentLoopTree::insert(SegmentLoop &arg_segmentloop, SegmentLoopTreeElement *arg_element)
{
	switch(arg_element->segmentloop.sequence.compare(arg_segmentloop.sequence))
	{
		case IS_SMALLER:// element is smaller than the inserted loop, means it has to go the the right
			if(arg_element->right != nullptr)
			{
				insert(arg_segmentloop, arg_element->right);
			}
			else
			{
				arg_element->right = new SegmentLoopTreeElement(arg_segmentloop);
			}
			break;
			
		case IS_LARGER:
			if(arg_element->left != nullptr)
			{
				insert(arg_segmentloop, arg_element->left);
			}
			else
			{
				arg_element->left = new SegmentLoopTreeElement(arg_segmentloop);
			}
			break;
#if DEBUG
		default://  inserting same element twice
			throw std::invalid_argument("SegmentLoopTree::insert() inserting same element twice.");
			break;
#endif //DEBUG
	}
}



/**
 * @brief Returns whether the tree is empty or not
 */
bool SegmentLoopTree::empty(void)
{
	return (this->root == nullptr);
}



/**
 * @brief Counts the elements in the subtree recursively
 */
size_t SegmentLoopTree::size(SegmentLoopTreeElement *arg_element)
{
	return (arg_element == nullptr) ? 0 : 1 + this->size(arg_element->left) + this->size(arg_element->right);
}



/**
 * @brief Counts the elements in the tree recursively
 */
size_t SegmentLoopTree::size(void)
{
	return this->empty() ? 0 : this->size(this->root);
}




/**
 * @brief Destructs the tree element recursively
 */
void SegmentLoopTree::clear(SegmentLoopTreeElement *arg_element)
{
	if(arg_element != nullptr)
	{
		this->clear(arg_element->left);
		this->clear(arg_element->right);
		
		delete arg_element;
	}
}



/**
 * @brief Destructs entire tree recursively
 */
void SegmentLoopTree::clear(void)
{
	this->clear(this->root);
}
