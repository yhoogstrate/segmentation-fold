/**
 * @file include/Pairing.hpp
 *
 * @date 2015-05-02
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



#ifndef PAIRINGTYPE_HPP
#define	PAIRINGTYPE_HPP


/**
 * @brief The canonical types of RNA-pairing
 * 
 * @date 2014-05-16
 */
enum PairingType : signed char
{
	AU = 0,
	CG = 1,
	GC = 2,
	UA = 3,
	GU = 4,
	UG = 5,
	None = -1
};

#endif	// PAIRINGTYPE_HPP
