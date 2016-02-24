/**
 * @file include/PairingPlus.hpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
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
 */


#ifndef PAIRINGPLUS_HPP
#define	PAIRINGPLUS_HPP


#include "main.hpp"

#include "Pair.hpp"
#include "PairingType.hpp"
#include "Position.hpp"
#include "Pairing.hpp"
#include "Sequence.hpp"



/**
 * @todo Merge Pair and Pairing into one object that has functionality of both by storing the two positions as iterators
 */
class PairingPlus: public Pairing
{
#if DEBUG
	private:
		void _check_order(void);// tests if position1 is indeed smaller than position2; t
#endif //DEBUG
		
	public:
		Position position1;
		Position position2;
		
		size_t size;/// The number of nucleotides in-between the pairing! [19,20] << has a size of 0
		
		PairingPlus(Position arg_position1, Position arg_position2);
		PairingPlus(Position arg_position1, Position arg_position2, size_t arg_nucleotides_inbetween);
		PairingPlus(Sequence &arg_sequence, Pair arg_pair);
		
		Nucleotide operator[](size_t);
};

#endif	// PAIRINGPLUS_HPP
