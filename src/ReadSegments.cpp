/**
 * @file src/ReadSegments.cpp
 *
 * @date 2015-08-06
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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
 * </PRE>
 */



#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctype.h>

#include <array>
#include <vector>

#include "Pair.hpp"
#include "Region.hpp"
#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"

#include "Direction.hpp"
#include "Sequence.hpp"
#include "Segment.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"

#include "ReadSegments.hpp"

#include <boost/foreach.hpp>
#include <boost/optional/optional.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>



using namespace boost;
using namespace boost::property_tree;


/**
 * @brief Initializes the ReadSegments class
 *
 * @date 2015-07-23
 */
ReadSegments::ReadSegments(std::string &arg_filename):
	filename(arg_filename)
{
	this->segments = nullptr;
}



/**
 * @brief Starts parsing of the <segments> section of the xml file
 *
 * @date 2015-07-23
 */
void ReadSegments::parse(SegmentTree &arg_segments)
{
	this->segments = (&arg_segments);
	
	this->parse(false);
}



/**
 * @brief Starts parsing of the <segments> and <rnas> sections of the xml file
 *
 * @date 2015-07-23
 */
void ReadSegments::parse(SegmentTree &arg_segments, std::vector<rna_example> &arg_examples)
{
	this->segments = (&arg_segments);
	this->rna_examples = (&arg_examples);
	
	this->parse(true);
}



/**
 * @brief Parses an XML file using boost xml library
 *
 * @date 2015-07-15
 */
void ReadSegments::parse(bool arg_parse_examples)
{
	std::ifstream ifs(this->filename.c_str(), std::ifstream::in);
	
	ptree xml_root;
	read_xml(ifs, xml_root);
	
	ptree xml_segments = xml_root.get_child("root.segments");
	this->parse_segments(xml_segments);
	
	if(arg_parse_examples)
	{
		ptree xml_examples = xml_root.get_child("root.rnas");
		this->parse_examples(xml_examples);
	}
	
	ifs.close();
}



/**
 * @brief Parses the <segments> section of the XML file
 *
 * @date 2015-07-23
 *
 * @section DESCRIPTION
 * <PRE>
 * For inversion the following needs to happen
 *    ------>
 * 5') ABCDE (3'    [p5 sequence]
 *     |  ||        [bonds]
 * 3') C  BA (5'    [p3 sequence]
 *    <------
 *
 * becomes:
 *
 *    ------>
 * 5') AB  C (3'    [p3 sequence]
 *     ||  |        [bonds]
 * 3') EDCBA (5'    [p5 sequence]
 *    <------
 *
 * Thus: p5, p3 and bonds are reversed and thereby p5 and p3 are swapped
 * for each other
 * </PRE>
 */
void ReadSegments::parse_segments(ptree &xml_segments)
{
	std::string str_true = "true";										///@todo make it static class member
	
	BOOST_FOREACH(ptree::value_type const & xml_segment, xml_segments)
	{
		if(xml_segment.first == "segment")
		{
			std::string id = xml_segment.second.get<std::string>("id");
			std::string p5 = xml_segment.second.get<std::string>("sequence_5prime");
			std::string bonds = xml_segment.second.get<std::string>("bonds");
			std::string p3 = xml_segment.second.get<std::string>("sequence_3prime");
			std::string energy = xml_segment.second.get<std::string>("energy");
			
			if(str_true.compare(xml_segment.second.get<std::string>("directions.five_prime")) == 0)
			{
				this->segments->insert(*(this->parse_segment(id, p5, bonds, p3, energy)));
			}
			if(str_true.compare(xml_segment.second.get<std::string>("directions.three_prime")) == 0)
			{
				// REVERSING THE SEQUENCE IS AN ILLUSION -- segments always go from 5' to 3', this doesn't change:
				// std::reverse(p5.begin(), p5.end());
				// std::reverse(bonds.begin(), bonds.end());
				// std::reverse(p3.begin(), p3.end());
				//
				// 5') xxxABCxxxh
				//     ||||||||| h
				// 3') xxxFEDxxxh
				//
				// We have:
				//    5') ABC (3'
				//    5') DEF (3'
				//
				// if we swap:
				//
				//    5') DEF (3'
				//    5') ABC (3'
				//
				// 5') xxxDEFxxxh
				//     ||||||||| h
				// 3') xxxCBAxxxh
				//
				
				std::reverse(p5.begin(), p5.end());
				std::reverse(bonds.begin(), bonds.end());
				std::reverse(p3.begin(), p3.end());
				
				this->segments->insert(*(this->parse_segment(id, p3, bonds, p5, energy)));// just swap their places
			}
		}
	}
}



/**
 * @brief Parses the <rna> section of the XML file
 *
 * @date 2015-07-15
 */
void ReadSegments::parse_examples(ptree &xml_examples)
{
	BOOST_FOREACH(ptree::value_type const & xml_rna, xml_examples)
	{
		if(xml_rna.first == "rna")
		{
			std::string title        = xml_rna.second.get<std::string>("title");
			std::string organism     = xml_rna.second.get<std::string>("organism", ""); //Default value is an empty string
			std::string sequence_str = xml_rna.second.get<std::string>("sequence");
			Sequence sequence        = Sequence(sequence_str);
			
			std::string dot_bracket  = "";
			
			BOOST_FOREACH(ptree::value_type const & xml_dot_bracket, xml_rna.second.get_child("structures.structure"))
			{
				if(xml_dot_bracket.first == "dot_bracket")
				{
					std::string tmp_dot_bracket_pattern = xml_dot_bracket.second.get<std::string>("<xmlattr>.type", "");
					if((tmp_dot_bracket_pattern == "segment") || (tmp_dot_bracket_pattern == "full" && dot_bracket == ""))
					{
						dot_bracket = xml_dot_bracket.second.data();
					}
				}
			}
			
			this->rna_examples->push_back(rna_example { title ,  organism ,  sequence , std::vector<Segment *>() , dot_bracket });
		}
	}
}



/**
 * @brief Parses a single segment (correct bonds) based on the XML data
 *
 * @date 2015-08-06
 */
Segment *ReadSegments::parse_segment(std::string arg_name, std::string arg_sequence_5p, std::string arg_bonds, std::string arg_sequence_3p, std::string arg_energy)
{
	std::vector<Pair> bonds = std::vector<Pair>();
	
	std::string abs_sequence_5p = arg_sequence_5p;
	std::string abs_sequence_3p = arg_sequence_3p;
	
	abs_sequence_5p.erase(std::remove(abs_sequence_5p.begin(), abs_sequence_5p.end(), ' '), abs_sequence_5p.end());
	abs_sequence_3p.erase(std::remove(abs_sequence_3p.begin(), abs_sequence_3p.end(), ' '), abs_sequence_3p.end());
	
	std::reverse(abs_sequence_3p.begin(), abs_sequence_3p.end());
	
	Sequence sequence_5p = Sequence(abs_sequence_5p);
	Sequence sequence_3p = Sequence(abs_sequence_3p);
	
	float energy = std::atof(arg_energy.c_str());
	
	int i = 1;
	int j = 1;
	
	unsigned int k;
	
	for(k = 0; k < arg_bonds.size(); k++)
	{
		if(arg_bonds[k] != ' ')
		{
			bonds.push_back(Pair {i, j});
			i = 0;
			j = 0;
		}
		
		if(arg_sequence_5p[k] != ' ')
		{
			i++;
		}
		if(arg_sequence_3p[k] != ' ')
		{
			j++;
		}
	}
	
	Segment *m = new Segment(arg_name, sequence_5p, bonds, sequence_3p, energy);
	this->segment_list.push_back(m);
	
	return m;
}



/**
 * @brief Removes all segments (memory safe)
 *
 * @date 2015-07-23
 */
void ReadSegments::clear(void)
{
	for(std::vector<Segment *>::iterator it = this->segment_list.begin(); it != this->segment_list.end(); ++it)
	{
		delete *it;///@note According to a valgrind test, brackets are not neccesairy
	}
	
	this->segment_list.clear();
}



/**
 * @brief Destructor
 *
 * @date 2015-07-23
 */
ReadSegments::~ReadSegments()
{
	this->clear();
}
