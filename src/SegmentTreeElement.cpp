/**
 * @file src/SegmentTreeElement.cpp
 *
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
#include "Segment.hpp"
#include "SegmentTreeElement.hpp"



/**
 * @brief Constructs a SegmentTreeElement by setting all pointers to nullptr and inserting the Segment object
 *
 */
SegmentTreeElement::SegmentTreeElement(Segment &arg_segment) :
	segment(arg_segment)
{
	this->link_5p_size_smaller = nullptr;
	this->link_5p_size_larger = nullptr;
	
	this->link_5p_segment_smaller = nullptr;
	this->link_5p_segment_larger = nullptr;
	
	
	this->link_3p_size_smaller = nullptr;
	this->link_3p_size_larger = nullptr;
	
	this->link_3p_segment_smaller = nullptr;
	this->link_3p_segment_larger = nullptr;
}



/**
 * @brief Destructor
 *
 */
SegmentTreeElement::~SegmentTreeElement()
{
	delete this->link_3p_segment_smaller;
	delete this->link_3p_segment_larger;
	
	delete this->link_3p_size_smaller;
	delete this->link_3p_size_larger;
	
	
	delete this->link_5p_segment_smaller;
	delete this->link_5p_segment_larger;
	
	delete this->link_5p_size_smaller;
	delete this->link_5p_size_larger;
}



/**
 * @brief Adds a Segment (as member of a SegmentTreeElement) to another element of the tree.
 *
 *
 * @todo use constants or static enum for cases 1, 2, 3 and 4
 */
void SegmentTreeElement::add_segment(Segment &arg_segment)
{
	char search_type = 1;
	
	this->add_segment(arg_segment, search_type);
}
void SegmentTreeElement::add_segment(Segment &arg_segment, char &arg_search_type)
{
	unsigned int i;
	
	Nucleotide self_n;
	Nucleotide arg_n;
	
	switch(arg_search_type)
	{
		case 1:
			if(this->segment.size(Direction::FivePrime) > arg_segment.size(Direction::FivePrime))
			{
				if(this->link_5p_size_smaller != nullptr)
				{
					this->link_5p_size_smaller->add_segment(arg_segment, arg_search_type);
				}
				else
				{
					SegmentTreeElement *element;
					element = new SegmentTreeElement(arg_segment);
					this->link_5p_size_smaller = element;
				}
			}
			else if(arg_segment.size(Direction::FivePrime) > this->segment.size(Direction::FivePrime))
			{
				if(this->link_5p_size_larger != nullptr)
				{
					this->link_5p_size_larger->add_segment(arg_segment, arg_search_type);
				}
				else
				{
					SegmentTreeElement *element;
					element = new SegmentTreeElement(arg_segment);
					this->link_5p_size_larger = element;
				}
			}
			else
			{
				this->add_segment(arg_segment, ++arg_search_type);
			}
			break;
			
		case 2:
			i = 0;
			
			while(i == 0 || (self_n == arg_n and i < this->segment.size(Direction::FivePrime)))
			{
				self_n = this->segment.get_nucleotide(Direction::FivePrime, i);
				arg_n = arg_segment.get_nucleotide(Direction::FivePrime, i);
				
				if(arg_n < self_n)
				{
					if(this->link_5p_segment_smaller != nullptr)
					{
						this->link_5p_segment_smaller->add_segment(arg_segment, arg_search_type);
					}
					else
					{
						SegmentTreeElement *element;
						element = new SegmentTreeElement(arg_segment);
						this->link_5p_segment_smaller = element;
					}
				}
				else if(arg_n > self_n)
				{
					if(this->link_5p_segment_larger != nullptr)
					{
						this->link_5p_segment_larger->add_segment(arg_segment, arg_search_type);
					}
					else
					{
						SegmentTreeElement *element;
						element = new SegmentTreeElement(arg_segment);
						this->link_5p_segment_larger = element;
					}
				}
				
				i++;
			}
			
			if(arg_n == self_n)
			{
				this->add_segment(arg_segment, ++arg_search_type);
			}
			
			break;
		case 3:
			if(this->segment.size(Direction::ThreePrime) > arg_segment.size(Direction::ThreePrime))
			{
				if(this->link_3p_size_smaller != nullptr)
				{
					this->link_3p_size_smaller->add_segment(arg_segment, arg_search_type);
				}
				else
				{
					SegmentTreeElement *element;
					element = new SegmentTreeElement(arg_segment);
					this->link_3p_size_smaller = element;
				}
			}
			else if(arg_segment.size(Direction::ThreePrime) > this->segment.size(Direction::ThreePrime))
			{
				if(this->link_3p_size_larger != nullptr)
				{
					this->link_3p_size_larger->add_segment(arg_segment, arg_search_type);
				}
				else
				{
					SegmentTreeElement *element;
					element = new SegmentTreeElement(arg_segment);
					this->link_3p_size_larger = element;
				}
			}
			else
			{
				this->add_segment(arg_segment, ++arg_search_type);
			}
			break;
		case 4:
			i = 0;
			
			while(i == 0 || (self_n == arg_n and i < this->segment.size(Direction::ThreePrime)))
			{
				self_n = this->segment.get_nucleotide(Direction::ThreePrime, i);
				arg_n = arg_segment.get_nucleotide(Direction::ThreePrime, i);
				
				if(arg_n < self_n)
				{
					if(this->link_3p_segment_smaller != nullptr)
					{
						this->link_3p_segment_smaller->add_segment(arg_segment, arg_search_type);
					}
					else
					{
						SegmentTreeElement *element;
						element = new SegmentTreeElement(arg_segment);
						this->link_3p_segment_smaller = element;
					}
				}
				else if(arg_n > self_n)
				{
					if(this->link_3p_segment_larger != nullptr)
					{
						this->link_3p_segment_larger->add_segment(arg_segment, arg_search_type);
					}
					else
					{
						SegmentTreeElement *element;
						element = new SegmentTreeElement(arg_segment);
						this->link_3p_segment_larger = element;
					}
				}
				
				i++;
			}
			
			if(arg_n == self_n)
			{
				std::string err = std::string("adding the following duplicate Segment multiple times to the SegmentTree:\n 5') " + arg_segment.get_sequence(Direction::FivePrime)->str() + "\n 3') " + arg_segment.get_sequence(Direction::ThreePrime)->str());
				
				throw std::invalid_argument(err.c_str());
			}
			break;
	}
}



/**
 * @brief Searches whether a segment exists in the tree
 *
 *
 * @todo use constants or static enum for cases 1, 2, 3 and 4
 */
Segment *SegmentTreeElement::search_segment(SubSequence &arg_sequence_5p, SubSequence &arg_sequence_3p)
{
	char search_type = 1;
	return this->search_segment(arg_sequence_5p, arg_sequence_3p, search_type);
}
Segment *SegmentTreeElement::search_segment(SubSequence &arg_sequence_5p, SubSequence &arg_sequence_3p, char &arg_search_type)
{
	Segment *m = nullptr;
	unsigned int i;
	
	Nucleotide self_n;
	Nucleotide arg_n;
	
	switch(arg_search_type)
	{
		case 1:///@todo use static constants for these cases
			if(this->segment.size(Direction::FivePrime) > arg_sequence_5p.size)
			{
				if(this->link_5p_size_smaller != nullptr)
				{
					m = this->link_5p_size_smaller->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
				}
			}
			else if(arg_sequence_5p.size > this->segment.size(Direction::FivePrime))
			{
				if(this->link_5p_size_larger != nullptr)
				{
					m = this->link_5p_size_larger->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
				}
			}
			else
			{
				m = this->search_segment(arg_sequence_5p, arg_sequence_3p, ++arg_search_type);
			}
			break;
			
		case 2:
			i = 0;
			
			/// @todo This loop can be improved!
			while(i == 0 || (arg_n == self_n && i < this->segment.size(Direction::FivePrime)))
			{
				self_n = this->segment.get_nucleotide(Direction::FivePrime, i);
				arg_n = arg_sequence_5p[i];
				
				if(arg_n < self_n)
				{
					if(this->link_5p_segment_smaller != nullptr)
					{
						m = this->link_5p_segment_smaller->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
					}
				}
				else if(arg_n > self_n)
				{
					if(this->link_5p_segment_larger != nullptr)
					{
						m = this->link_5p_segment_larger->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
					}
				}
				
				i++;
			}
			
			if(arg_n == self_n)
			{
				m = this->search_segment(arg_sequence_5p, arg_sequence_3p, ++arg_search_type);
			}
			
			break;
		case 3:
			if(this->segment.size(Direction::ThreePrime) > arg_sequence_3p.size)
			{
				if(this->link_3p_size_smaller != nullptr)
				{
					m = this->link_3p_size_smaller->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
				}
			}
			else if(arg_sequence_3p.size > this->segment.size(Direction::ThreePrime))
			{
				if(this->link_3p_size_larger != nullptr)
				{
					m = this->link_3p_size_larger->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
				}
			}
			else
			{
				m = this->search_segment(arg_sequence_5p, arg_sequence_3p, ++arg_search_type);
			}
			break;
		case 4:
			i = 0;
			
			while(i == 0 || (arg_n == self_n && i < this->segment.size(Direction::ThreePrime)))
			{
				self_n = this->segment.get_nucleotide(Direction::ThreePrime, i);
				arg_n = arg_sequence_3p[i];
				
				if(arg_n < self_n)
				{
					if(this->link_3p_segment_smaller != nullptr)
					{
						m = this->link_3p_segment_smaller->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
					}
				}
				else if(arg_n > self_n)
				{
					if(this->link_3p_segment_larger != nullptr)
					{
						m = this->link_3p_segment_larger->search_segment(arg_sequence_5p, arg_sequence_3p, arg_search_type);
					}
				}
				
				i++;
			}
			
			if(arg_n == self_n)
			{
				m = &this->segment;
			}
			
			break;
	}
	
	return m;
}



/**
 * @brief Counts recursively the number of elements in the tree
 *
 */
size_t SegmentTreeElement::size(void)
{
	size_t n = 1;														// Adds itself to the total
	
	if(this->link_5p_size_smaller != nullptr)
	{
		n += this->link_5p_size_smaller->size();
	}
	if(this->link_5p_size_larger != nullptr)
	{
		n += this->link_5p_size_larger->size();
	}
	
	if(this->link_5p_segment_smaller != nullptr)
	{
		n += this->link_5p_segment_smaller->size();
	}
	if(this->link_5p_segment_larger != nullptr)
	{
		n += this->link_5p_segment_larger->size();
	}
	
	if(this->link_3p_size_smaller != nullptr)
	{
		n += this->link_3p_size_smaller->size();
	}
	if(this->link_3p_size_larger != nullptr)
	{
		n += this->link_3p_size_larger->size();
	}
	
	if(this->link_3p_segment_smaller != nullptr)
	{
		n += this->link_3p_segment_smaller->size();
	}
	if(this->link_3p_segment_larger != nullptr)
	{
		n += this->link_3p_segment_larger->size();
	}
	
	return n;
}
