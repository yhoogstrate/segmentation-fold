/**
 * @file src/Utils/file_exists.cpp
 * @date 22-apr-2015
 * @author Youri Hoogstrate
 * @section DESCRIPTION
 * Stores the file_exists() function.
 *
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
 * @brief Checks if a file exists.
 * @date 12-mar-2013
 * @param filename The filename of which presency should be checked.
 * @return true if the file exists, false if it doesn't exist.
 */
bool file_exists(const char *filename)
{
	std::ifstream ifile(filename);
	return ifile;
}
