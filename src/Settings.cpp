/**
 * @file src/Settings.cpp
 *
 * @date 2015-07-20
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



#include "Utils/utils.hpp"

#include "Pair.hpp"
#include "Region.hpp"
#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"

#include "ReadSegments.hpp"
#include "ScoringTree.hpp"
#include "ReadData.hpp"

#include "Settings.hpp"



/**
 * @brief Construct of the Settings class; default parameters should be initialized in this function.
 *
 * @date 2015-06-22
 *
 * @param local_argc Local argc (number of commandline parameters) parameter which should be passed through to the global argc
 * @param local_argv Local argv (commandline parameters) parameter which should be passed through to the global argv
 * @param arg_sequence An empty Sequence object that will be set by this class
 */
Settings::Settings(int arg_argc, char **arg_argv, Sequence &arg_sequence) :
	obj_sequence(arg_sequence),
	argc(arg_argc),
	argv(arg_argv)
{
	this->minimal_hairpin_length = 3;
	this->segment_prediction_functionality = true;
	this->num_threads = 0;
	
	this->segment_filename = std::string();
	
	this->parse_arguments();
	this->get_segments_file();
}



/**
 * @brief Prints usage
 *
 * @date 2015-07-23
 */
void Settings::print_usage(void)
{
	fprintf(stderr, "Usage: " PACKAGE_NAME " -s [SEQUENCE]\n");
	fprintf(stderr, "Usage: " PACKAGE_NAME " -f [FASTA_FILE]\n");
	fprintf(stderr, "   * Note: If FASTA_FILE and SEQUENCE are not provided,\n           the program will read from STDIN.\n\n\n");
	fprintf(stderr, "The following parameters can be used:\n");
	fprintf(stderr, "  -s SEQUENCE                Specific RNA SEQUENCE (overrules -f)\n");
	fprintf(stderr, "  -f FASTA_FILE              Location to FASTA_FILE that contains the sequence\n\n");
	fprintf(stderr, "  -p                  [1/0]  Enable/disable segment prediction functionality\n\n");
	fprintf(stderr, "  -h HAIRPINSIZE      [1,N}  Minimum hairpin size, default: 3\n");
	fprintf(stderr, "  -x SEGMENTS_XML_FILE       Use custom  \"segments.xml\"-syntaxed file\n\n");
	fprintf(stderr, "  -t NUM_THREADS      [0,N}  Number of threads; 0 = maximum available, default: 3.\n\n");
	fprintf(stderr, "  -V                         Shows the version and license\n");
	fprintf(stderr, "\n\n");
	fprintf(stderr, "If you encounter problems with this software, please send bug-reports to:\n   <" PACKAGE_BUGREPORT ">\n\n");
	
	this->run_print_usage = true;
}



/**
 * @brief Prints the package version, git commit sha1, build-type and license
 *
 * @date 2015-04-23
 */
void Settings::print_version(void)
{
#if DEBUG
#define BUILD_TYPE_STRING " (debug)"
#else
#define BUILD_TYPE_STRING " (release)"
#endif

	printf("[Version]\n  " PACKAGE_STRING GIT_SHA1_STRING BUILD_TYPE_STRING "\n\n");
	printf("[License]\n  GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n\n  This is free software: you are free to change and redistribute it.\n  There is NO WARRANTY, to the extent permitted by law.\n\n  Written by Youri Hoogstrate.\n");
	
	this->run_print_version = true;
}



/**
 * @brief Parses the commandline parameters.
 *
 * @date 2015-12-01
 */
void Settings::parse_arguments(void)
{
	// Ensure argument parsing always starts from the 1st element
	// This allows parsing arguments multiple timesin one program
	optind = 1;
	
	bool proceed = true;
	char c;
	size_t i;
	
	if(this->argc > 1 && strcmp(this->argv[1], "--version") == 0)
	{
		this->argv[1] = (char *) "-V\0";
	}
	
	while((c = getopt(this->argc, this->argv, "+h:f:s:p:x:t:V")) > 0 && proceed)
	{
		switch(c)
		{
			case 'h':													// option -h for minimum hairpin length
				for(i = 0; i < strlen(optarg); i++)///@todo Validate whether we are converting a true integer
				{
					if(!isdigit(optarg[i]))
					{
						this->print_usage();
					}
				}
				sscanf(optarg, "%d", &this->minimal_hairpin_length);	// TODO use atoi?
				break;
			case 'f':
				FILE *stream;
				stream = fopen(optarg, "r");
				
				if(stream != nullptr)
				{
					this->parse_sequence_from_stream(stream);
					fclose(stream);
				}
				else
				{
					throw std::invalid_argument("Can't open file \"" + std::string(optarg) + "\".");
				}
				break;
			case 's':
				this->obj_sequence = Sequence(optarg);
				break;
			case 'p':
				this->segment_prediction_functionality = (optarg[0] == '1');
				break;
			case 'x':
				if(file_exists(optarg))
				{
					this->segment_filename = std::string(optarg);
				}
				else
				{
					throw std::invalid_argument("Can't open file \"" + std::string(optarg) + "\".");
				}
				break;
			case 't':
				for(i = 0; i < strlen(optarg); i++)
				{
					if(!isdigit(optarg[i]))
					{
						this->print_usage();
					}
				}
				sscanf(optarg, "%d", &this->num_threads);
				break;
			case 'V':
				this->print_version();
				proceed = false;
				break;
			default:
				this->print_usage();
				proceed = false;
				break;
		}
	}
	
	if(this->obj_sequence.empty() && this->run_print_version == false && this->run_print_usage == false)
	{
		printf("Please insert your RNA sequence:\n");
		this->parse_sequence_from_stream(stdin);
	}
}



/**
 * @brief Parses the first FASTA line from a filestream (or stdin)
 *
 * @date 2015-06-22
 *
 * @param stream A stream for reading, like a file handle or stdin.
 *
 * @todo find a MAXLEN - probably related to max size of int
 */
void Settings::parse_sequence_from_stream(FILE *stream)
{
	char buffer;
	
	while(((buffer = fgetc(stream)) != EOF))// and (obj_sequence.size() < MAXLEN - 2))
	{
		if(buffer == '>')												// comment line for FASTA format start with '>'
		{
			if(this->obj_sequence.empty())
			{
				do
				{
					buffer = fgetc(stream);
				}
				while(buffer != '\n');
			}
			else
			{
				break;
			}
		}
		else
		{
			if(buffer == 'A' || buffer == 'U' || buffer == 'T' || buffer == 'C' || buffer == 'G' ||
					buffer == 'a' || buffer == 'u' || buffer == 't' || buffer == 'c' || buffer == 'g')
			{
				this->obj_sequence.push_back(buffer);
			}
		}
	}
}



/**
 * @brief Finds directories that should contain the segmentation-fold share data.
 *
 * @date 2015-06-23
 *
 * @return The directories found in environment variable $XDG_DATA_DIRS.
 */
std::vector<std::string> Settings::get_share_directories(void)
{
	unsigned int i;
	
	std::vector<std::string> directories = std::vector<std::string>();
	directories.push_back(DATA_DIR "/share/");
	
	char *env_datadirs;
	env_datadirs = getenv("XDG_DATA_DIRS");
	
	char *env_datadirs_buf;
	if(env_datadirs != nullptr)
	{
		env_datadirs_buf = NULL;
		
		for(i = 0; i < strlen(env_datadirs); i++)
		{
			if(env_datadirs_buf == nullptr)
			{
				env_datadirs_buf = &env_datadirs[i];
			}
			
			if(env_datadirs[i] == ':' or env_datadirs[i] == '\0')
			{
				env_datadirs[i] = '\0';
				directories.push_back(std::string(env_datadirs_buf));
				env_datadirs_buf = nullptr;
			}
		}
	}
	
	directories.push_back("../share/");
	directories.push_back("share/");
	directories.push_back(SEGMENTS_PATH "/");
	directories.push_back("~/.local/share/");
	
	return directories;
}



/**
 * @brief Find the segments.xml file
 *
 * @section DESCRIPTION
 * The file is usually found in "/usr/local/share/[package_name]
 *
 * @date 2013-07-28
 */
void Settings::get_segments_file(void)
{
	// Check whether it was not already being set by an argument
	if(this->segment_filename.empty())									// If the filename was not already argumented
	{
		std::vector<std::string> data_directories = this->get_share_directories();
		std::string file_name = std::string();
		
		bool search_file = true;
		
		for(
			std::vector<std::string>::iterator it = data_directories.begin();
			it != data_directories.end() && search_file;
			++it
		)
		{
			file_name = *it + PACKAGE_NAME "/" SEGMENTS_FILE;
			
			if(file_exists(file_name.c_str()))
			{
				this->segment_filename = file_name;
				search_file = false;
			}
			else
			{
				file_name = *it + SEGMENTS_FILE;
				
				if(file_exists(file_name.c_str()))
				{
					this->segment_filename = file_name;
					search_file = false;
				}
			}
		}
		
		// if it still hasn't been found
		if(this->segment_filename.empty())
		{
			throw std::invalid_argument("Couldn't load the \"" SEGMENTS_FILE "\" file");
		}
	}
}
