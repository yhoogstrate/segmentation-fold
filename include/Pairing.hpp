/**
 * @file include/Pairing.hpp
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



#ifndef PAIRING_HPP
#define	PAIRING_HPP


#include "PairingType.hpp"


/**
 * @brief Container of two Nucleotides
 */
class Pairing
{
	public:
		PairingType type;
		
		Pairing(Nucleotide arg_n1, Nucleotide arg_n2);
		
		void init(Nucleotide arg_n1, Nucleotide arg_n2);
		
		bool is_canonical(void);
};


#endif	// PAIRING_HPP
