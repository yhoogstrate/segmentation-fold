/**
 * @date 2015-04-22
 *
 * @file include/ScoringTree.hpp
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

#ifndef SCORINGTREE_HPP
#define	SCORINGTREE_HPP


#include "ScoringTreeElement.hpp"


template <class T_key, class T_value>
class ScoringTree
{
	private:
		ScoringTreeElement<T_key, T_value> *root;
		void insert(T_key arg_key, T_value &arg_value, ScoringTreeElement<T_key, T_value> *arg_parent);
		T_value *search(T_key &arg_key, ScoringTreeElement<T_key, T_value> *arg_parent);
		void clear(ScoringTreeElement<T_key, T_value> *arg_parent);
		
	public:
		ScoringTree();
		~ScoringTree();
		
		void insert(T_key arg_key, T_value &arg_value);
		T_value *search(T_key &arg_key);
		void clear();
};


#endif	// SCORINGTREE_HPP
