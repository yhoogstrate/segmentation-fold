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
a bioinformatics application that predicts RNA 2D-structure
with an extended version of the Zuker algorithm. This modification contains a new "structure element"
named a *segment* is capable of folding a pre-defined substructure with **multiple** canonical
or non-canonical pairings.

This allows folding of more complex structures like the K-turns, which
are also part of the implemented free energy tables. These thermodynamic
parameters (free Gibbs energy levels) have been estimated using a
computational approach and lack therefore accuracy.

# INSTALLATION #

## Prerequisites ##

Segmentation-fold does **not** install on most systems of itself because
it depends on two additional libraries and an installation library.

	cmake
	boost library (-test)
	boost library (-xml)

In Ubuntu and Debian you can install these packages with the following command:

	sudo apt-get install libboost-test-dev libboost-filesystem-dev

In Arch linux you can install these packages with the folllowing command:

	sudo pacman -S git cmake boost boost-libs

## Recommended packages ##

To create the corresponding documentation, you should have installed the
following package:

	doxygen (>= 1.8.3)

The doxygen package is version specific because of the Markdown
support implemented in 1.8.3 and above.

To automatically make the synax of the code tidy you should have
installed the package:

	astyle

You can install the recommended packages in Ubuntu or Debian with:

	sudo apt-get install doxygen astyle

You can install the recommanded packages in Arch with:

	sudo pacman -S doxygen astyle

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
between debug and release are that (1) debug contains certain checks on
variables and (2) has no optimizated (but fast) compilation. Therefore,
debug will result in a slower executing, but more secure binary.

## Get the documentation ##

Run the following commands to get the latest version of the
documentation:

	cmake .
	make readme
	make doc

The first command, "cmake ." creates the doxygen description file from
its template and sets e.g. the package version variable. Then
"make readme" creates the latests version of README.md based on the Markdown
files present in this project. Then "make doc" creates the documentation
from all the doxygen comments within the code, including the latest
README.md as main page.

# USAGE #

	Usage: segmentation-fold -h minimum_hairpin_length -s [SEQUENCE]
	Usage: segmentation-fold -h minimum_hairpin_length -f [FASTA_FILE]
	   * Note: If FASTA_FILE and SEQUENCE are not provided,
	           the program will read from STDIN.

	The following parameters can be used:
	  -h HAIRPINSIZE             Minimum hairpin size
	  -f FASTA_FILE              Location to FASTA_FILE
	  -s SEQUENCE                Specific RNA SEQUENCE (overrules -f)
	  -k [1/0]                   Enable/disable K-turn predictions
	  -m MOTIFS_FILE             Use custom  "motifs.xml"-syntaxed file
	  -t [0/1]                   Use novel (=0) or original (=1) traceback
	  -v                         Shows the version and license

## galaxy ##

segmentation-fold will become available as galaxy tool.
# AUTHORS #

All people that have contributed to the development of segmentation-fold are:


# CONTRIBUTING #

We encourage users, researchers and programmers to contribute to this free and open source project. This can be achieved by reporting bugs and commiting code to the github repository https://github.com/yhoogstrate/segmentation-fold. To streamline and archive communication in an univocal way, we encourge contributors to only use this channel, github, to contribute to segmentation-fold.

If any developer wants to contribute to segmentation-fold, please make use the following documentation standard: http://www.stack.nl/~dimitri/doxygen/index.html

All source-code is formatted using Asteric Style: http://astyle.sourceforge.net/ The corresponding configuration file is available under 'share/.astylerc'. If you want to contribute to the code you should have it installed to 'beautify' the code, by simply running:

	make tidy

Cmake will make use of the correct syntax file and beautify all C++ and hpp files in src/, include/ and test/.
