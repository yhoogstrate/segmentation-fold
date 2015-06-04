/**
 * @file src/ScoringTree.cpp
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



#include "main.hpp"
#include "Utils/utils.hpp"

#include <array>

#include "Pair.hpp"
#include "Direction.hpp"
#include "Nucleotide.hpp"
#include "Sequence.hpp"
#include "Segment.hpp"



#include "ScoringTree.hpp"



/**
 * @brief
 * 
 * @date 2015-05-02
 */
template <class T_key, class T_value>
ScoringTree<T_key, T_value>::ScoringTree()
{
	this->root = nullptr;
}



/**
 * @brief
 * 
 * @date 2015-05-02
 */
template <class T_key, class T_value>
ScoringTree<T_key, T_value>::~ScoringTree()
{
	this->clear();
}



/**
 * @brief
 * 
 * @date 2015-05-02
 */
template <class T_key, class T_value>
void ScoringTree<T_key, T_value>::insert(T_key arg_key, T_value &arg_value)
{
	if(this->root != nullptr)
	{
		this->insert(arg_key, arg_value, this->root);
	}
	else
	{
		this->root = new ScoringTreeElement<T_key, T_value> {arg_key, arg_value, nullptr, nullptr};
	}
}



/**
 * @brief
 * 
 * @date 2015-05-02
 */
template <class T_key, class T_value>
void ScoringTree<T_key, T_value>::insert(T_key arg_key, T_value &arg_value, ScoringTreeElement<T_key, T_value> *arg_parent)
{
	if(arg_key < arg_parent->key)
	{
		if(arg_parent->smaller != nullptr)
		{
			this->insert(arg_key, arg_value, arg_parent->smaller);
		}
		else
		{
			arg_parent->smaller = new ScoringTreeElement<T_key, T_value> {arg_key , arg_value, nullptr, nullptr};
		}
	}
	else if(arg_key > arg_parent->key)
	{
		if(arg_parent->larger != nullptr)
		{
			this->insert(arg_key, arg_value, arg_parent->larger);
		}
		else
		{
			arg_parent->larger = new ScoringTreeElement<T_key, T_value> {arg_key , arg_value, nullptr, nullptr};
		}
	}
#if DEBUG
	else
	{
		// The algorithm should never insert a segment ort a score at the same Pair twice
		throw std::invalid_argument("ScoringTree::ScoringTree(): Trying to overwrite a ScringTree element.\n");
	}
#endif //DEBUG
}



/**
 * @brief
 * 
 * @date 2015-05-02
 */
template <class T_key, class T_value>
T_value *ScoringTree<T_key, T_value>::search(T_key &arg_key)
{
	return this->search(arg_key, this->root);
}



template <class T_key, class T_value>
T_value *ScoringTree<T_key, T_value>::search(T_key &arg_key, ScoringTreeElement<T_key, T_value> *arg_parent)
{
	if(arg_parent != nullptr)
	{
		if(arg_key == arg_parent->key)
		{
			return &arg_parent->value;
		}
		else if(arg_key < arg_parent->key)
		{
			return this->search(arg_key, arg_parent->smaller);
		}
		else
		{
			return this->search(arg_key, arg_parent->larger);
		}
	}
	else
	{
		return nullptr;
	}
}



template <class T_key, class T_value>
void ScoringTree<T_key, T_value>::clear()
{
	this->clear(this->root);
}



template <class T_key, class T_value>
void ScoringTree<T_key, T_value>::clear(ScoringTreeElement<T_key, T_value> *arg_parent)
{
	if(arg_parent != nullptr)
	{
		this->clear(arg_parent->smaller);
		this->clear(arg_parent->larger);
		
		delete arg_parent;
	}
}



// Explicit instantiation:
template class ScoringTree<Pair, unsigned int>;
template class ScoringTree<Pair, Segment>;
