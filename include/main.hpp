/**
 * @file include/main.hpp
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

#include <limits.h>
#include <float.h>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

// http://stackoverflow.com/a/7266266 << template way of doing it
// http://stackoverflow.com/a/28338942
//#define N_INFINITY  FLT_MAX
#define MAX_OF(type) \
    (((type)(~0LLU) > (type)((1LLU<<((sizeof(type)<<3)-1))-1LLU)) ? (long long unsigned int)(type)(~0LLU) : (long long unsigned int)(type)((1LLU<<((sizeof(type)<<3)-1))-1LLU))
#define MIN_OF(type) \
    (((type)(1LLU<<((sizeof(type)<<3)-1)) < (type)1) ? (long long int)((~0LLU)-((1LLU<<((sizeof(type)<<3)-1))-1LLU)) : 0LL)

#define N_INFINITY         (float) MAX_OF(float)

#define BOUND              -1
#define UNBOUND            -2
#define NOT_YET_CALCULATED -3

#define MISCLOOP_PRELOG                    0

#define MISCLOOP_ASYMETRIC_INTENRAL_LOOP   1

#define MISCLOOP_FM_ARRAY_0                2
#define MISCLOOP_FM_ARRAY_1                3
#define MISCLOOP_FM_ARRAY_2                4
#define MISCLOOP_FM_ARRAY_3                5

#define MISCLOOP_MULTIBRANCH_OFFSET        6
#define MISCLOOP_MULTIBRANCH_FB_PENALTY    7
#define MISCLOOP_MULTIBRANCH_HELIX_PENALTY 8

#define MISCLOOP_EFN2_OFFSET               9
#define MISCLOOP_EFN2_FB_PENALTY           10
#define MISCLOOP_EFN2_HELIX_PENALTY        11

#define MISCLOOP_AU_PENALTY                12
#define MISCLOOP_GGG_U_PENALTY             13
#define MISCLOOP_C_HAIRPIN_SLOPE           14
#define MISCLOOP_C_HAIRPIN_INTERCEPT       15
#define MISCLOOP_C_HAIRPIN_OF_3            16

#define MISCLOOP_INTERMOLECULAR_INIT       17
#define MISCLOOP_GAIL_RULE                 18

#define DOTBRACKET__PAIRING_LEFT  '('
#define DOTBRACKET__PAIRING_RIGHT ')'
#define DOTBRACKET__NO_PAIRING    '.'

#define IS_EQUAL                  0
#define IS_SMALLER                1
#define IS_LARGER                 2


#endif	// MAIN_HPP
