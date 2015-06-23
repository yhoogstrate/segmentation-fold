/**
 * @file include/Settings.hpp
 *
 * @date 2014-04-17
 *
 * @author Youri Hoogstrate
 * @author Lisa Yu
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


#ifndef SETTINGS_HPP
#define	SETTINGS_HPP


#include "main.hpp"
#include <config.hpp>

/**
 * @brief Takes care of the input, arguments and stores the settings; actual parsing of files within arguments is not part of this class.
 *
 * @section DESCRIPTION
 * Handling the arguments, settings and input.
 *
 * @date 2015-06-22
 */
class Settings
{
	private:
		///@todo reconsider whether to store the variables rather than argument them directly - will they ever be reused?
		int argc;
		char **argv;
		
		Sequence &obj_sequence;
		
		void parse_arguments(void);
		void parse_sequence_from_stream(FILE *stream);
		
		std::vector<std::string> get_share_directories(void);
		void get_segments_file(void);
		
		// Print functions
		void print_usage(void);
		void print_version(void);
	public:
		Settings(int arg_argc, char **arg_argv, Sequence &arg_sequence);
		
		int minimal_hairpin_length;
		bool segment_prediction_functionality;
		std::string segment_filename;
		
		bool run_print_version = false;
		bool run_print_usage = false;
};


#endif	// SETTINGS_HPP
