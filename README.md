[![Build Status](https://travis-ci.org/yhoogstrate/segmentation-fold.svg?branch=master)](https://travis-ci.org/yhoogstrate/segmentation-fold)
___

	 _______  _______  _______  _______  _______  _       _________ _______ __________________ _______  _       
	(  ____ \(  ____ \(  ____ \(       )(  ____ \( (    /|\__   __/(  ___  )\__   __/\__   __/(  ___  )( (    /|
	| (    \/| (    \/| (    \/| () () || (    \/|  \  ( |   ) (   | (   ) |   ) (      ) (   | (   ) ||  \  ( |
	| (_____ | (__    | |      | || || || (__    |   \ | |   | |   | (___) |   | |      | |   | |   | ||   \ | |
	(_____  )|  __)   | | ____ | |(_)| ||  __)   | (\ \) |   | |   |  ___  |   | |      | |   | |   | || (\ \) |
	      ) || (      | | \_  )| |   | || (      | | \   |   | |   | (   ) |   | |      | |   | |   | || | \   |
	/\____) || (____/\| (___) || )   ( || (____/\| )  \  |   | |   | )   ( |   | |   ___) (___| (___) || )  \  |
	\_______)(_______/(_______)|/     \|(_______/|/    )_)   )_(   |/     \|   )_(   \_______/(_______)|/    )_)
	                                                                                                            
	                                                     _____                                                  
	                                                    (_____)                                                 
	                                                                                                            
	                                       _______  _______  _        ______                                    
	                                      (  ____ \(  ___  )( \      (  __  \                                   
	                                      | (    \/| (   ) || (      | (  \  )                                  
	                                      | (__    | |   | || |      | |   ) |                                  
	                                      |  __)   | |   | || |      | |   | |                                  
	                                      | (      | |   | || |      | |   ) |                                  
	                                      | )      | (___) || (____/\| (__/  )                                  
	                                      |/       (_______)(_______/(______/                                   
___

# INTRODUCTION #

[segmentation-fold](https://github.com/yhoogstrate/segmentation-fold) is
a bioinformatics application that predicts RNA 2D-structure with an
extended version of the Zuker algorithm. This modification contains a
new "structure element" named a *segment* and is capable of folding a
pre-defined substructure with **multiple** canonical or non-canonical
pairings.

This allows folding of more complex structures like the K-turns, which
are also part of the implemented free energy tables. These thermodynamic
parameters (free Gibbs energy levels) have been estimated using a
computational approach and therefore lack accuracy.

# INSTALLATION #

## Prerequisites ##

Segmentation-fold does **not** install on most systems of itself because
it depends on two additional libraries and an installation library.

	cmake
	boost library (-system)
	boost library (-test)
	boost library (-xml)

In Ubuntu and Debian you can install these packages with the following command:

	sudo apt-get install cmake libboost-system-dev libboost-test-dev libboost-filesystem-dev

In Arch linux you can install these packages with the folllowing command:

	sudo pacman -S git cmake boost boost-libs

## Recommended packages ##

To create the corresponding documentation, you should have installed the
following package:

	doxygen (>= 1.8.3)

The doxygen package is version specific because Markdown support is
implemented in 1.8.3 and above.

To automatically make (new) code styled similar to the programming style
using in this package, you should have installed the package:

	astyle

You can install the recommended packages in Ubuntu or Debian with:

	sudo apt-get install doxygen astyle

You can install the recommanded packages in Arch with:

	sudo pacman -S doxygen astyle

Details on running doxygen and astyle are given in the sections
[Get the documentation](https://github.com/yhoogstrate/segmentation-fold#get-the-documentation)
and
[CONTRIBUTING](https://github.com/yhoogstrate/segmentation-fold#contributing).

## Build and install ##

After having the prerequisites installed, you can downadload, compile
and install segmentation-fold with following commands:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release .
	 $ make
	 $ make test
	 $ sudo make install

In case you do not have administrator rights on a particular machine,
you can install as follows:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=~/.local .
	 $ make
	 $ make test
	 $ make install

You will then find the binary in the directory:
	<home directory>/.local/bin/

If you want to run segmentation-fold for anything else than development,
you should compile with the -DCMAKE_BUILD_TYPE=release flag because by
default compiling is done with the 'debug' flag. The two differences
between debug and release are that (1) debug contains certain
(unneccesairy if all code behaves as expected) checks on variables and
(2) has no optimized, but faster, compilation. Therefore, debug mode
will produce a binary that has a slower performance, but a higher level
of safety.

If you are using a compiled version of segmentation-fold and you don't
known the original compilation settings, you should run:

	segmentation-fold --version

This will report the version number, the git commit-SHA1 and the build
type.

## Get the documentation ##

Run the following commands to get the latest version of the
documentation:

	cmake .
	make readme
	make doc

The first command, "cmake ." creates the doxygen description file from
its template and sets e.g. the package version variable. Then "make
readme" creates the latests version of README.md based on the Markdown
files present in this project. Then "make doc" creates the documentation
from all the doxygen comments within the code, including the latest
README.md as main page.

# USAGE #

	Usage: segmentation-fold -s [SEQUENCE]
	Usage: segmentation-fold -f [FASTA_FILE]
	   * Note: If FASTA_FILE and SEQUENCE are not provided,
	           the program will read from STDIN.
	
	
	The following parameters can be used:
	  -s SEQUENCE                Specific RNA SEQUENCE (overrules -f)
	  -f FASTA_FILE              Location to FASTA_FILE that contains the sequence
	
	  -p                  [1/0]  Enable/disable segment prediction functionality
	
	  -h HAIRPINSIZE      [1,N}  Minimum hairpin size, default: 3
	  -x SEGMENTS_XML_FILE       Use custom  "segments.xml"-syntaxed file
	
	  -t NUM_THREADS      [0,N}  Number of threads; 0 = maximum available, default: 3.
	
	  -V                         Shows the version and license
	
	
	If you encounter problems with this software, please send bug-reports to:
	   <https://github.com/yhoogstrate/segmentation-fold/issues>

## galaxy ##

segmentation-fold will become available for galaxy. The wrapper will be
available at:
[https://github.com/ErasmusMC-Bioinformatics/segmentation_fold_galaxy_wrapper](https://github.com/ErasmusMC-Bioinformatics/segmentation_fold_galaxy_wrapper)
# AUTHORS #

All people that have contributed to the development of segmentation-fold are:


# CONTRIBUTING #

We encourage users, researchers and programmers to contribute to this
free and open source project. This can be achieved by reporting bugs and
commiting code to the github repository
https://github.com/yhoogstrate/segmentation-fold. To streamline and
archive communication in an univocal way, we encourge contributors to
only use this channel, Github, to contribute to segmentation-fold.
To contribute to segmentation-fold, please make use the following
documentation standard:
[http://www.stack.nl/~dimitri/doxygen/index.html](http://www.stack.nl/~dimitri/doxygen/index.html).

## Tidy ##

All source-code is formatted using Asteric Style: [http://astyle.sourceforge.net/](http://astyle.sourceforge.net/).
The corresponding configuration file is available at '[share/.astylerc](https://raw.githubusercontent.com/yhoogstrate/segmentation-fold/master/share/.astylerc)'.
If you want to contribute to the code you should have it installed to 'beautify' the code, by simply running:

	make tidy

Cmake will make use of the correct syntax file and beautify all C++ and hpp files in
[src/](https://github.com/yhoogstrate/segmentation-fold/tree/master/src),
[include/](https://github.com/yhoogstrate/segmentation-fold/tree/master/include),
and
[test/](https://github.com/yhoogstrate/segmentation-fold/tree/master/test).

## Tests ##

To simplify the reviewing process of submitted code the project contains
unit and functional tests. These tests have to be passed in order to get
a positive review. These tests also inspect memory leaks using valgrind.
To run the test on (your copy of the) code before doing a pull request, run:

	cmake .
	make clean
	make readme
	make
	make check
	ctest -V -T memcheck

This will re-build the readme, re-compile the code, and does testing with and
without memory leak checking. If you can't get it working but you still believe
your change is worth submitting, don't worry. Whenever you do a pull request
TravisCI will automatically run the tests for you.
