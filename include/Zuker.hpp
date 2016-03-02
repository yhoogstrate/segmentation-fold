/**
 * @file include/Zuker.hpp
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * <PRE>
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2016 Youri Hoogstrate
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
 * </PRE>
 */



#ifndef ZUKER_HPP
#define	ZUKER_HPP


#include "main.hpp"
#include "Utils/utils.hpp"


#define V_MATRIX 1
#define W_MATRIX 2
#define WM_MATRIX 3



#include "ZukerTraceback.hpp"
#include "GibbsFreeEnergy.hpp"


/**
 * @brief Implements the Zuker algorithm with additional segment-folding features.
 *
 * @section DESCRIPTION
 * This file contains the Zuker class.
 * The program implements the Zuker's minimum Gibbs free energy energy
 * model with Segments/K-turns functionality for RNA secondary structure
 * prediction.
 */
class Zuker: public GibbsFreeEnergy
{
		friend class Test_Zuker;
		
	private:
		Settings &settings;
		
		std::vector<traceback_jump> traceback_stack;
		int traceback_stacktop;
		unsigned int folded_segments;
		
		Position sequence_begin;
		
	public:
		Zuker(Settings &arg_settings, Sequence &arg_sequence, ReadData &arg_thermodynamics);
		
		DotBracket dot_bracket;
		
		// Energy per structure functions:
		float energy(void);// Fills the V and W matrices
		
		// Energy functions:
		float v(Pair &p1, PairingPlus &p1p);
		float w(Pair &p1);
		float wm(Pair &p1, PairingPlus &p1p);
		
		// Trace-back related:
		void traceback(void);
		void traceback_push(traceback_jump arg_jump);
		bool traceback_pop(unsigned int &i, unsigned int &j, char &matrix);
		
		// Output functions
		void print_2D_structure(void);
		
		// Energy matrices
		ScoringMatrix<float> vij;//paired matrix
		ScoringMatrix<float> wij;//unpaired matrix
		ScoringMatrix<float> wmij;//multiloop matrix
		
		// Traceback matrices
		ScoringMatrix<traceback_jump> tij_v;
		ScoringMatrix<traceback_jump> tij_w;
		ScoringMatrix<traceback_jump> tij_wm;
		ScoringMatrix<SegmentTraceback *> sij;
		
#if DEBUG
		// Functions only useful for plotting and debugging
		void _print_sij(unsigned int matrix_length);
		void _print_vij(unsigned int matrix_length);
		void _print_wij(unsigned int matrix_length);
#endif //DEBUG
};



///@brief Friend class of Zuker that allows testing its private members
class Test_Zuker: public Zuker
{
	public:
		using Zuker::vij;
		using Zuker::Zuker;
};



#endif	// ZUKER_HPP
