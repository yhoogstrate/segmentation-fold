/**
 * @file src/SubSequence.cpp
 *
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



#include "main.hpp"

#include "Nucleotide.hpp"
#include "Position.hpp"
#include "SubSequence.hpp"



/**
 * @brief Stores positions and sets PairingType
 *
 */
SubSequence::SubSequence(Position arg_position1, Position arg_position2):
	position1(arg_position1),
	position2(arg_position2),
	size((size_t)(arg_position2 - arg_position1  + 1))
{
	this->init();
}



/**
 * @brief Stores positions and sets PairingType - If size is known apriori, it is useless to re-calculate it and thus save it immediately
 *
 */
SubSequence::SubSequence(Position arg_position1, Position arg_position2, size_t arg_size):
	position1(arg_position1),
	position2(arg_position2),
	size(arg_size)
{
	this->init();
}


/**
 * @brief Sets the PairingType
 *
 */
void SubSequence::init()
{
#if DEBUG
	this->_check_order();
#endif //DEBUG
}


/**
 * @brief
 *
 *
 * @todo Decide whether or not to return a Position or a Nucleotide
 */
Nucleotide SubSequence::operator[](size_t arg_position)
{
#if DEBUG
	if(arg_position >= this->size)
	{
		throw std::invalid_argument("SubSequence::operator[]: Out of bound SubSequence");
	}
#endif //DEBUG
	return *(this->position1 + (long) arg_position);
}



#if DEBUG
/**
 * @brief checks whether the order of the position is correct and throws an exception otherwise
 */
void SubSequence::_check_order(void)
{
	if(this->position2 < this->position1)
	{
		throw std::invalid_argument("Position 2 of SubSequence is smaller than Position 1");
	}
}
#endif //DEBUG
