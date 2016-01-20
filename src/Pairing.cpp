/**
 * @file src/Pairing.cpp
 *
 * @date 2016-01-20
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

#include "Nucleotide.hpp"
#include "Pairing.hpp"

/**
 * @date 2014-03-10
 */
Pairing::Pairing(Nucleotide arg_n1, Nucleotide arg_n2)
{
	this->init(arg_n1, arg_n2);
}



/**
 * @date 2016-01-20
 */
void Pairing::init(Nucleotide arg_n1, Nucleotide arg_n2)
{
	switch(arg_n1)
	{
		case Nucleotide::A:
			switch(arg_n2)
			{
				case Nucleotide::U:
					this->type = PairingType::AU;
					return void();
					break;
				default:														// Makes the compiler happy..
					break;
			}
			break;
			
		case Nucleotide::C:
			switch(arg_n2)
			{
				case Nucleotide::G:
					this->type = PairingType::CG;
					return void();
					break;
				default:														// Makes the compiler happy..
					break;
			}
			break;
			
		case Nucleotide::G:
			switch(arg_n2)
			{
				case Nucleotide::U:
					this->type = PairingType::GU;
					return void();
					break;
				case Nucleotide::C:
					this->type = PairingType::GC;
					return void();
					break;
				default:														// Makes the compiler happy..
					break;
			}
			break;
			
		case Nucleotide::U:
			switch(arg_n2)
			{
				case Nucleotide::A:
					this->type = PairingType::UA;
					return void();
					break;
				case Nucleotide::G:
					this->type = PairingType::UG;
					return void();
				default:														// Makes the compiler happy..
					break;
					break;
			}
			break;
	}
	
	this->type = PairingType::None;
}



/**
 * @brief Return whether the pairing is a canonical RNA pairing.
 *
 * @date 2014-03-10
 */
bool Pairing::is_canonical(void)
{
	return (this->type != PairingType::None);
}
