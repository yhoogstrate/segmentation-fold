/**
 * @file include/ReadSegments.hpp
 *
 * @date 2015-12-07
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



#include <boost/property_tree/ptree.hpp>

using namespace boost::property_tree;



/**
 * @brief Structure for an (example) RNA that contains a segment
 *
 * @date 2015-04-22
 */
struct rna_example
{
	std::string title;
	std::string organism;
	Sequence sequence;
	
	std::vector<Segment *> segments;
	//std::vector<SegmentLoop *> segmentloops;
	std::string dot_bracket_pattern;
};



/**
 * @brief Parses the segment XML file
 *
 * @date 2015-12-07
 */
class ReadSegments
{
	private:
		std::string &filename;
		
		SegmentTree *segments;
		std::vector<Segment *> segment_list;
		std::vector<SegmentLoop *> segmentloop_list;
		std::vector<rna_example> *rna_examples;
		
		void parse(bool arg_parse_examples);
		
		void parse_segments(ptree &xml_motifs);
		void parse_examples(ptree &xml_examples);
		
		Segment *parse_segment(std::string arg_name, std::string arg_sequence_5p, std::string arg_bonds, std::string arg_sequence_3p, std::string arg_energy);
		
	public:
		ReadSegments(std::string &arg_filename);
		
		void parse(SegmentTree &arg_segments, SegmentLoopTree &arg_segmentloops);
		void parse(SegmentTree &arg_segments, SegmentLoopTree &arg_segmentloops, std::vector<rna_example> &arg_examples);
		
		void clear(void);
		~ReadSegments();
};
