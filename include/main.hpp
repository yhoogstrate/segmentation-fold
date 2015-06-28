/**
 * @file include/main.hpp
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

#ifndef MAIN_HPP
#define	MAIN_HPP


#include <config.hpp>

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdarg.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <string>
#include <vector>


#define BOUND              -1
#define UNBOUND            -2
#define NOT_YET_CALCULATED -3

#define DOTBRACKET__PAIRING_LEFT  '('
#define DOTBRACKET__PAIRING_RIGHT ')'
#define DOTBRACKET__NO_PAIRING    '.'

#endif	// MAIN_HPP
