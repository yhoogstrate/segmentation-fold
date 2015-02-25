/**
 * @file include/Direction.hpp
 *
 * @date 25-feb-2015
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2014 Youri Hoogstrate
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



#ifndef DIRECTION_HPP
#define	DIRECTION_HPP


/**
 * @brief Describes a direction (5' to 3' or vice versa) in at most one byte
 * 
 * @date 13-mar-2014
 */
enum struct Direction : bool
{
	/*
	Backwards , ThreePrime = true,
	Forwards  , FivePrime  = false
	*/
	
	FivePrime  = true,													// 5' => 3'
	ThreePrime = false													// 3' => 5'
};


#endif																	// DIRECTION_HPP

