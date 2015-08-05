/**
 * @file src/Segment.cpp
 *
 * @date 2015-05-02
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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
 * </PRE>
 */



#include "main.hpp"

#include <array>

#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "Segment.hpp"



/**
 * @brief Constructor of the Segment class.
 *
 * @section DESCRIPTION
 * For the following motif:
 *
 * 5') ...CCCCC....
 *        | | |    .
 * 3') ...A U G....
 *
 * Alignment:   Bonds:
 * C-A          0-2
 * C
 * C-U          2-1
 * C
 * C-G          4-0
 *
 * We denote:
 * name="arbitrary name"
 * sequence_5p="CCCCC"
 * sequence_3p="GUA"   <- this inversion is correct, because then it's annotated from 5' to 3'
 * bonds={{0,2},{2,1},{4,0}}
 *
 * @param arg_name Name of the segment
 * @param arg_sequence_5p The segments sequence located closer to the 5' end of the RNA
 * @param arg_sequence_3p The segments sequence located closer to the 3' end of the RNA
 * @param arg_bonds A description of the bonds between the the two sequences, as integers starting from 0
 *
 * @date 2014-03-21
 */
Segment::Segment(std::string arg_name, Sequence arg_sequence_5p, std::vector <Pair> arg_bonds, Sequence arg_sequence_3p, float arg_gibbs_free_energy):
	name(arg_name),
	sequence_5p(arg_sequence_5p),
	sequence_3p(arg_sequence_3p),
	bonds(arg_bonds),
	gibbs_free_energy(arg_gibbs_free_energy)
{
	this->reset_traceback();
}



/**
 * @brief Resets the traceback to the first bond of the segment
 *
 * @date 2015-04-16
 */
void Segment::reset_traceback(void)
{
	this->it = this->bonds.rbegin();
}



/**
 * @brief Return the size of the 5' or 3' sequence
 *
 * @section DESCRIPTION
 * In the following example:
 *
 * 5') ...CCCCC...  <- Direction::FivePrime sequence
 *        | | |
 * 3') ...A AAA.... <- Direction::ThreePrime sequence
 *
 * Direction::FirePrime will return 5 (CCCCC = 5 nucleotides)
 * Direction::TreePrime will return 4 (AAAA = 4 nucleotides)
 *
 * @date 2015-05-02
 * */
size_t Segment::size(Direction direction)
{
	return (direction == Direction::FivePrime) ? this->sequence_5p.size() : this->sequence_3p.size();
}
size_t Segment::size(Direction &direction)
{
	return (direction == Direction::FivePrime) ? this->sequence_5p.size() : this->sequence_3p.size();
}



/**
 * @brief Finds the next offset for the traceback algorithm
 *
 * @section DESCRIPTION
 * For the following examplesequence
 *
 * 5') ...CCCCC...  <- Direction::FivePrime sequence
 *        | | |
 * 3') ...A U G.... <- Direction::ThreePrime sequence
 *
 * We have the following corresponding information:
 * Alignment:   Bonds:      Traceback:
 * C-A          0-2         -1,-1
 * C
 * C-U          2-1         -3,-2
 * C
 * C-G          4-0         -5,-3
 *
 * Where bonds is what's located within this segment class.
 * The traceback list is the information that should be returned when
 * the corresponding value in Bonds is found.
 *
 * @param i Reference to variable i that should be set
 * @param j Reference to variable i that should be set
 * @return Whether the return was VALID - if false is returned, the traceback is being reset to it's origin and will return true again.
 *
 * @date 2014-04-15
 *
 * @todo Find a way to be more efficient with signed/unsigned
 */
bool Segment::pop(signed int &i, signed int &j)
{
	if(this->it == this->bonds.rend())
	{
		this->reset_traceback();
		
		return false;
	}
	else
	{
		Pair pair = *this->it;
		
		i = - (this->sequence_5p.size() - pair.first);
		j = -                             pair.second - 1;
		
		this->it++;///@todo see if this iterator can move 5 lines up and run tests whether it doesn't break
		
		return true;
	}
}



/**
 * @brief Returns the free energy that corresponds to the segment
 *
 * @date 2015-08-03
 *
 * @todo check if it is possible to do by sending the reference and see whether that's faster
 */
const float Segment::get_gibbs_free_energy(void)
{
	return this->gibbs_free_energy;
}



/**
 * @brief Returns a particular nucleotide in one of the 2 sequences of the segment at a given direciton
 *
 * @param direction Direction::FivePrime is used if the sequence is located closest to the 5' end of the RNA sequence; Direction::ThreePrime otherwise
 * @param i Location (starting from 0, ending at n-1) in the sequence of interest
 *
 * @date 2015-04-23
 */
Nucleotide Segment::get_nucleotide(Direction direction, unsigned int &i)
{
	return (direction == Direction::FivePrime) ? this->sequence_5p[i] : this->sequence_3p[i];
}



/**
 * @brief Returns one of the 2 sequences of the segment at a given direciton
 *
 * @param direction Direction::FivePrime is used if the sequence is located closest to the 5' end of the RNA sequence; Direction::ThreePrime otherwise
 *
 * @date 2015-04-16
 */
Sequence *Segment::get_sequence(Direction direction)
{
	return (direction == Direction::FivePrime) ? &this->sequence_5p : &this->sequence_3p;
}
