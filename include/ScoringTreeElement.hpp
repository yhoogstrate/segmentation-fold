/**
 * @date 2015-06-22
 *
 * @file include/ScoringTreeElement.hpp
 *
 * @todo Figure out whether std::map may not be more appropriate
 * A discussion on this can be found at
 * <http://stackoverflow.com/questions/5288320/why-is-stdmap-implemented-as-red-black-tree>
 * suggesting that instead of the std::map's red black tree the AVL tree
 * may also be a good choice.
 *
 * @section DESCRIPTION
 * The scoring matrices are always filled at every position. The segment
 * objects are only located at certain positions; where the segments are
 * indeed found in the sequence. Therefore a sparsers tree structure
 * enhances performance.
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

#ifndef SCORINGTREEELEMENT_HPP
#define	SCORINGTREEELEMENT_HPP

/**
 * @brief Tree element for ScoringTree
 * 
 * @date 2014-04-20
 */
template<class T_key, class T_value>
struct ScoringTreeElement
{
	T_key key;
	T_value &value;
	
	ScoringTreeElement<T_key, T_value> *smaller;
	ScoringTreeElement<T_key, T_value> *larger;
};


#endif	// SCORINGTREEELEMENT_HPP


// Explicit instantiation:
template struct ScoringTreeElement<Pair, unsigned int>;
template struct ScoringTreeElement<Pair, Segment>;
