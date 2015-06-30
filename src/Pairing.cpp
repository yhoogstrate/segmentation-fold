/**
 * @file src/Pairing.cpp
 *
 * @date 2014-03-10
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
 *
 * @author Youri Hoogstrate
 */
Pairing::Pairing(Nucleotide &arg_n1, Nucleotide &arg_n2)
{
	this->n1 = &arg_n1;
	this->n2 = &arg_n2;
	
	this->init(arg_n1, arg_n2);
}



/**
 * @date 04-apr-2015
 *
 * @todo check whether a multiplication and a switch can be more efficient
 */
void Pairing::init(Nucleotide &arg_n1, Nucleotide &arg_n2)
{
	if(arg_n1 == Nucleotide::A && arg_n2 == Nucleotide::U)
	{
		this->type = PairingType::AU;
	}
	else if(arg_n1 == Nucleotide::C && arg_n2 == Nucleotide::G)
	{
		this->type = PairingType::CG;
	}
	else if(arg_n1 == Nucleotide::G && arg_n2 == Nucleotide::C)
	{
		this->type = PairingType::GC;
	}
	else if(arg_n1 == Nucleotide::U && arg_n2 == Nucleotide::A)
	{
		this->type = PairingType::UA;
	}
	else if(arg_n1 == Nucleotide::G && arg_n2 == Nucleotide::U)
	{
		this->type = PairingType::GU;
	}
	else if(arg_n1 == Nucleotide::U && arg_n2 == Nucleotide::G)
	{
		this->type = PairingType::UG;
	}
	else
	{
		this->type = PairingType::None;
	}
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
