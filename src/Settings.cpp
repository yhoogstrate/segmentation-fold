/**
 * @file src/Settings.cpp
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



#include <libgen.h>
#include <limits.h>

#include "Utils/utils.hpp"

#include "Pair.hpp"
#include "Region.hpp"
#include "Nucleotide.hpp"
#include "Pairing.hpp"
#include "PairingPlus.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"
#include "SegmentLoop.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"
#include "SegmentLoopTree.hpp"

#include "ReadSegments.hpp"
#include "ReadData.hpp"

#include "Settings.hpp"



/**
 * @brief Construct of the Settings class; default parameters should be initialized in this function.
 *
 *
 * @param local_argc Local argc (number of commandline parameters) parameter which should be passed through to the global argc
 * @param local_argv Local argv (commandline parameters) parameter which should be passed through to the global argv
 * @param arg_sequence An empty Sequence object that will be set by this class
 */
Settings::Settings(int arg_argc, char **arg_argv, Sequence &arg_sequence) :
	argc(arg_argc),
	argv(arg_argv),
	obj_sequence(arg_sequence)
{
	this->minimal_hairpin_length = 3;
	this->segment_prediction_functionality = true;
	this->num_threads = 0;
	
	this->segment_filename = std::string();
	
	this->proceed_with_folding = true;
	
	this->parse_arguments();
}



/**
 * @brief Prints usage
 */
void Settings::print_usage(bool as_error)
{
	this->proceed_with_folding = false;
	
	std::ostream &stream = as_error ? std::cerr : std::cout;
	
	stream << "Usage: " PACKAGE_NAME " -s [SEQUENCE]\n";
	stream << "       " PACKAGE_NAME " -f [FASTA_FILE]\n";
	stream << "   * Note: If FASTA_FILE and SEQUENCE are not provided,\n";
	stream << "           the program will read from STDIN.\n";
	stream << "\n\n";
	stream << "The following parameters can be used:\n";
	stream << "  -s SEQUENCE       Specific RNA SEQUENCE (overrules -f)\n";
	stream << "  -f FASTA_FILE     Path of FASTA_FILE containing sequence(s)\n";
	stream << "  -p                Enable/disable segment functionality           [1/0]\n";
	stream << "  -H HAIRPINSIZE    Minimum hairpin size, default: 3               [1,N}\n";
	stream << "  -x SEGMENTS_XML   Use custom  \"segments.xml\"-syntaxed file\n";
	stream << "  -t NUM_THREADS    Number of threads; 0 = maximum available,      [0,N}\n";
	stream << "                    default: 3 \n\n";
	stream << "  -h, --help        Display this help and exit\n";
	stream << "  -V, --version     Show version and license\n";
	stream << "  -X, --default-xml Show path to default \"segments.xml\" on\n";
	stream << "                    system\n";
	stream << "\n\n";
	stream << "If you encounter problems with this software, please report it at:\n";
	stream << "   <" PACKAGE_BUGREPORT ">\n\n";
}



/**
 * @brief Prints the package version, git commit sha1, build-type and license
 */
void Settings::print_version(void)
{
	this->proceed_with_folding = false;
	
#if DEBUG
#define BUILD_TYPE_STRING " (debug)"
#else
#define BUILD_TYPE_STRING " (release)"
#endif
	
	std::cout << "[Version]\n  " PACKAGE_STRING GIT_SHA1_STRING BUILD_TYPE_STRING "\n\n";
	std::cout << "[License]\n  GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n\n";
	std::cout << "  This is free software: you are free to change and redistribute it.\n";
	std::cout << "  There is NO WARRANTY, to the extent permitted by law.\n\n";
	std::cout << "  Copyright (C) 2012-2017  Youri Hoogstrate.\n";
}


void Settings::print_default_xml(void)
{
	this->proceed_with_folding = false;
	
	/* convenient for development
		printf("scanning directories:\n");
		std::vector<std::string> data_directories = this->get_share_directories();
		for(std::vector<std::string>::iterator it = data_directories.begin(); it != data_directories.end(); ++it)
		{
			std::cout << "\t" << *it << PACKAGE_NAME << "/\n";
		}
		std::cout << "\n";
	*/
	
	this->get_segments_file();
	char *real_path = realpath(this->segment_filename.c_str(), nullptr);
	std::cout << real_path << "\n";
	
	free(real_path);
}



/**
 * @brief Parses the commandline parameters.
 */
void Settings::parse_arguments(void)
{
	// Ensure argument parsing always starts from the 1st element
	// This allows parsing arguments multiple times during one program (e.g. functional testing)
	optind = 1;
	
	bool proceed_parsing_arguments = true;//proceed with parsing arguments
	
	int c;
	size_t i;
	
	if(this->argc > 1)
	{
		if(strcmp(this->argv[1], "--version") == 0)
		{
			this->argv[1] = (char *) "-V\0";
		}
		else if(strcmp(this->argv[1], "--default-xml") == 0)
		{
			this->argv[1] = (char *) "-X\0";
		}
		else if(strcmp(this->argv[1], "--help") == 0)
		{
			this->argv[1] = (char *) "-h\0";
		}
	}
	
	// 'So to distinguish them, getopt provides a mechanism. All the options that require argument will be preceded by a : (colon).'
	while((c = getopt(this->argc, this->argv, "+H:f:s:p:x:t:hVX")) > 0 && proceed_parsing_arguments)
	{
		switch(c)
		{
			case 'H':							// option -H for minimum hairpin length
				for(i = 0; i < strlen(optarg); i++)
				{
					if(!isdigit(optarg[i]))
					{
						proceed_parsing_arguments = false;
						break;
					}
				}
				
				if(proceed_parsing_arguments == false)
				{
					this->print_usage(true);
					throw std::invalid_argument("Invalid argument (-" + std::string(1, (char) c) + "): " + std::string(optarg));
				}
				else
				{
					sscanf(optarg, "%d", &this->minimal_hairpin_length);	// TODO use atoi?
				}
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
					this->print_usage(true);
					throw std::invalid_argument("Invalid argument (-" + std::string(1, (char) c) + "): can't open file \"" + std::string(optarg) + "\"");
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
					this->print_usage(true);
					throw std::invalid_argument("Invalid argument (-" + std::string(1, (char) c) + "): can't open file \"" + std::string(optarg) + "\"");
				}
				break;
			case 't':
				for(i = 0; i < strlen(optarg); i++)
				{
					if(!isdigit(optarg[i]))
					{
						proceed_parsing_arguments = false;
						break;
					}
				}
				
				if(proceed_parsing_arguments == false)
				{
					this->print_usage(true);
					throw std::invalid_argument("Invalid argument (-" + std::string(1, (char) c) + "): " + std::string(optarg));
				}
				else
				{
					sscanf(optarg, "%d", &this->num_threads);
				}
				break;
			case 'V':
				proceed_parsing_arguments = false;
				this->print_version();
				break;
			case 'X':
				proceed_parsing_arguments = false;
				this->print_default_xml();
				break;
			default:
				proceed_parsing_arguments = false;
				this->print_usage(false);
				break;
		}
	}
	
	if(this->proceed_with_folding)
	{
		this->get_segments_file();
		
		if(this->obj_sequence.empty())
		{
			printf("Please insert your RNA sequence:\n");
			this->parse_sequence_from_stream(stdin);
		}
	}
}



/**
 * @brief Parses the first FASTA line from a filestream (or stdin)
 *
 * @param stream A stream for reading, like a file handle or stdin.
 *
 * @todo find a MAXLEN - probably related to max size of int
 */
void Settings::parse_sequence_from_stream(FILE *stream)
{
	int buffer;
	
	while(((buffer = fgetc(stream)) != EOF))
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
				this->obj_sequence.push_back((char) buffer);
			}
		}
	}
}



/**
 * @brief Finds directories that should contain the segmentation-fold share data.
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
		env_datadirs_buf = nullptr;
		
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
	
	
	directories.push_back(this->readlink_self() + "/../share/");
	
	return directories;
}


/**
 * @brief Returns path of executable in a mem safe way as std::string
 *
 * @static
 */
std::string Settings::readlink_self()
{
	char buffer[PATH_MAX + 1];
	ssize_t len = ::readlink("/proc/self/exe", buffer, PATH_MAX);
	if(len != -1)
	{
		buffer[len] = '\0';
		return std::string(buffer);
	}
	else
	{
		return std::string("");
	}
}


/**
 * @brief Find the segments.xml file
 *
 * @section DESCRIPTION
 * The file is usually found in "/usr/local/share/[package_name]
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
