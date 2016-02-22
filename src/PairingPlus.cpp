/**
 * @file src/PairingPlus.cpp
 *
 * @date 2016-01-21
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
#include "PairingPlus.hpp"



/**
 * @brief Stores positions and sets PairingType
 */
PairingPlus::PairingPlus(Position arg_position1, Position arg_position2):
	PairingPlus(arg_position1, arg_position2, (size_t)(arg_position2 - arg_position1 - 1))
{
}



/**
 * @brief Stores positions and sets PairingType - If size is known apriori, it is useless to re-calculate it and thus save it immediately
 */
PairingPlus::PairingPlus(Position arg_position1, Position arg_position2, size_t arg_size):
	Pairing(*arg_position1 , *arg_position2),
	position1(arg_position1),
	position2(arg_position2),
	size(arg_size)
{
#if DEBUG
	if(this->position2 <= this->position1)
	{
		throw std::invalid_argument("Position 2 of Pairing is smaller than Position 1");
	}
#endif //DEBUG
}
