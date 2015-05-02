/**
 * @file include/Segment.hpp
 * 
 * @date 2015-05-02
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


#ifndef SEGMENT_HPP
#define SEGMENT_HPP


#include <array>


/**
 * @brief This object represents a *non-canonical RNA pairing segment* which requires to be surrounded by canonical pairings, like the K-turn.
 * @date 2015-04-20
 * @section DESCRIPTION
 *
 * The following can be converted into a dot-bracket format:
 * ACUUG
 * | | |
 * A U G
 *
 * ACUUG***AUG
 * ( ( (   )))
 *
 * The following can not be moddeled in a dot bracket; the C and the U of the upper sequence are bound the the U.
 *
 * ACUUG
 * |\| |
 * A U G
 *
 * ACUUG***AUG
 * (?( (   )))
 *
 * Using a derived bpseq format you could still model this:
 *
 *  1  A 11
 *  2  C 10
 *  3  U 10
 *  4  U 0
 *  5  G 9
 *  6  * 0
 *  7  * 0
 *  8  * 0
 *  9  G 5
 *  10 U 3 2   <-- order should here be important
 *  11 A 1
 *
 *
 * Consequently, the input of the bonds format should be in such a way that it may include multiple bonds per nucleotide, in a human readable format.
 * Therefore, we propose to upload segments in the same format:
 *
 *   1 A -1
 *   2 C -2
 *   3 U -2
 *   4 U  0
 *   5 G -3
 *
 *  -3 G  5
 *  -2 U  3  2   <-- order should here be important
 *  -1 A  1
 *
 * and mirrored:
 *
 * A U G
 * |\| |
 * ACUUG
 *
 *   1 A -1 -3
 *   2 C -3
 *   3 U -5
 *  -5 U  3
 *  -4 G  0
 *  -3 A  1
 *  -2 U  0
 *  -1 G  1
 *
 * @date 2014-04-15
 */
class Segment
{
	private:
		std::vector<Pair>::reverse_iterator it;
		
	public:
		Segment(std::string arg_name, Sequence arg_segment5p, std::vector<Pair> arg_bonds, Sequence arg_segment3p, float arg_gibbs_free_energy);
		
		std::string name;
		Sequence sequence_5p;
		Sequence sequence_3p;
		std::vector<Pair> bonds;
		
		bool pop(int &i, int &j);
		void reset_traceback(void);
		
		float gibbs_free_energy;
		float get_gibbs_free_energy(void);
		
		Nucleotide get_nucleotide(Direction direction, unsigned int &i);// @todo operator: segment[Direction][arg_i]
		Sequence  *get_sequence(Direction direction);					// @todo operator: segment[Direction]
		
		size_t size(Direction &direction);
		size_t size(Direction direction);
};

#endif	// SEGMENT_HPP
