/**
 * @file src/Utils/min.cpp
 * @date 11-sept-2013
 * @author Youri Hoogstrate
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2015 Youri Hoogstrate
 *
 * This file is part of segmentation-fold and originally taken from
 * yh-kt-fold.
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



#include "Utils/utils.hpp"



/**
 * @brief Finds the minimal value between x and y
 * @date 05-nov-2012
 * @param x any integer
 * @param y any integer
 * @return the minimum of both the values x and y
 * @todo Use macro or inline instead?
 */
int min(int x, int y)
{
	return (x < y ? x : y);
}



/**
 * @brief Finds the minimal value between x and y
 * @date 30-dec-2014
 * @param x any integer
 * @param y any integer
 * @return the minimum of both the values x and y
 * @todo Use macro or inline instead?
 *//*
float min(float x, float y)
{
	return (x < y ? x : y);
}*/
