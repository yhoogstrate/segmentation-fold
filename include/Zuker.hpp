/**
 * @file include/Zuker.hpp
 *
 * @date 2015-07-13
 *
 * @author Youri Hoogstrate
 * @author Lisa Yu
 *
 * @section LICENSE
 * <PRE>
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
 * </PRE>
 */



#ifndef ZUKER_HPP
#define	ZUKER_HPP


#include "main.hpp"
#include "Utils/utils.hpp"


/**
 * @brief Jumping element for the traceback function
 *
 * @date 2013-10-02
 */
struct traceback_jump
{
	int i;
	int j;
	bool jump_to_v_path;
};

struct traceback_jump_pair
{
	Pair pair;
	bool jump_to_v_path;
};



#include "GibbsFreeEnergy.hpp"


/**
 * @brief Implements the Zuker algorithm with additional segment-folding features.
 *
 * @section DESCRIPTION
 * This file contains the Zuker class.
 * The program implements the Zuker's minimum Gibbs free energy energy
 * model with Segments/K-turns functionality for RNA secondary structure
 * prediction.
 *
 * @date 2015-12-01
 */
class Zuker: public GibbsFreeEnergy
{
	private:
		Settings &settings;
		
		std::vector<traceback_jump> traceback_stack;
		int traceback_stacktop;
		
		Position sequence_begin;
		
	public:
		Zuker(Settings &arg_settings, Sequence &arg_sequence, ReadData &arg_thermodynamics);
		
		DotBracket dot_bracket;
		
		// Energy per structure functions:
		float energy(void);// Fills the V and W matrices
		
		// Energy functions:
		float v(Pair &p1, PairingPlus &p1p);
		float w(Pair &p1);
		
		inline float energy_bifurcation(Region &region);
		
		// Trace-back related:
		void traceback(void);
		void traceback_push(int i, int j, bool pick_from_v_path);
		bool traceback_pop(int *i, int *j, bool *pick_from_v_path);
		
		// Output functions
		void print_2D_structure(void);
		
		ScoringMatrix<int> pij;// Pathmatrix
		ScoringMatrix<int> qij;// Pathmatrix <if from bifurcation>
		
		ScoringMatrix<float> vij;
		ScoringMatrix<float> wij;
		
		ScoringMatrix<char> pathmatrix_corrected_from;					//@todo benchmark performance of std::vector<std::vector<bool>> pathmatrix_corrected_from; or make argument between fast and sparse
		ScoringMatrix<Pair> loopmatrix;									//@todo benchmark performance of std::vector<std::vector<bool>> pathmatrix_corrected_from; or make argument between fast and sparse
		
		ScoringMatrix<Segment *> nij2;									//@todo benchmark performance using std::map<Pair &>
		
#if DEBUG
		// Functions only useful for plotting and debugging
		void _print_pathmatrix_corrected_from(unsigned int matrix_length);
		void _print_loopmatrix(unsigned int matrix_length);
		void _print_nij(unsigned int matrix_length);
		void _print_pij(unsigned int matrix_length);
		void _print_vij(unsigned int matrix_length);
		void _print_wij(unsigned int matrix_length);
#endif //DEBUG
};

#endif	// ZUKER_HPP
