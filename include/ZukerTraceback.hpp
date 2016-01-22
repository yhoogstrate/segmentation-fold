/**
 * @file include/ZukerTraceback.hpp
 *
 * @date 2016-01-20
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
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
 * </PRE>
 */



#ifndef ZUKERTRACEBACK_HPP
#define	ZUKERTRACEBACK_HPP


#include "main.hpp"
#include "Pair.hpp"


///@date 2016-01-20
struct traceback_jump2
{
	bool store_pair;// is True, position (i,j) have to be stored as pair: i:'(', j:')'
	
	// target.first == target.second implies bifurcation between (i, target.first) & (target.first + 1, j)
	// target.first < 0, end of line; stop traceback -> preferably if target.first > n (sequence size)
	Pair target;
};

#endif	// ZUKERTRACEBACK_HPP
