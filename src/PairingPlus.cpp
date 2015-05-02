/**
 * @file src/PairingPlus.cpp
 *
 * @date 2015-04-29
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



#include "main.hpp"

#include "Nucleotide.hpp"
#include "Position.hpp"
#include "PairingPlus.hpp"



/**
 * @brief Stores positions and sets PairingType
 *
 * @date 2015-04-29
 *
 * @todo check what access by reference does/will do?
 */
PairingPlus::PairingPlus(Position arg_position1, Position arg_position2):
	size(arg_position2 - arg_position1 + 1)
{
	this->position1 = arg_position1;
	this->position2 = arg_position2;
	
	this->init();
}


/**
 * @brief Sets the PairingType
 * 
 * @date 2015-04-29
 */
void PairingPlus::init()
{
	if(*this->position1 == Nucleotide::A && *this->position2 == Nucleotide::U)
	{
		this->type = PairingType::AU;
	}
	else if(*this->position1 == Nucleotide::C && *this->position2 == Nucleotide::G)
	{
		this->type = PairingType::CG;
	}
	else if(*this->position1 == Nucleotide::G && *this->position2 == Nucleotide::C)
	{
		this->type = PairingType::GC;
	}
	else if(*this->position1 == Nucleotide::U && *this->position2 == Nucleotide::A)
	{
		this->type = PairingType::UA;
	}
	else if(*this->position1 == Nucleotide::G && *this->position2 == Nucleotide::U)
	{
		this->type = PairingType::GU;
	}
	else if(*this->position1 == Nucleotide::U && *this->position2 == Nucleotide::G)
	{
		this->type = PairingType::UG;
	}
	else
	{
		this->type = PairingType::None;
	}
}

/**
 * @brief Returns whether the pairing is a canonical RNA pairing.
 *
 * @date 2015-04-29
 */
bool PairingPlus::is_canonical(void)
{
	return (this->type != PairingType::None);
}


///@todo create a class SubSequence(PairingPlus) and make this member
Nucleotide PairingPlus::operator[](size_t arg_position)
{
	return *(this->position1 + arg_position);
}

