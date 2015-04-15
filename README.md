# segmentation-fold #

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

<https://github.com/yhoogstrate/segmentation-fold/>

Segmentation-fold is a bioinformatics application that predicts RNA 2D-structure
with an extended version of the Zuker algorithm. This modification contains a new "structure element"
named a *segment* is capable of folding a pre-defined substructure with **multiple** canonical
or non-canonical pairings.

	This allows folding of more complex structures like the K-turns, which are also part of
	the implemented free energy tables. These thermodynamic
parameters (free Gibbs energy levels) have been estimated using a
computational approach and lack therefore accuracy.


### INSTALLATION ###

#### Prerequisites ####

Segmentation-fold does **not** install on most systems of itself because
it depends on two additional libraries and an installation library.

	cmake
	boost library (-test)
	boost library (-xml)

In ubuntu you can install them very easily with the following command:

	sudo apt-get install libboost-test-dev libboost-filesystem-dev

#### Build and install ####

Please follow the instructions given in the '[INSTALL](https://github.com/yhoogstrate/segmentation-fold/blob/master/INSTALL)' file, e.g. by
using one of the following commands:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release .
	 $ make
	 $ make test
	 $ sudo make install

In case you do not have administrator rights on a particular machine, you can
install as follows:

	 $ git clone https://github.com/yhoogstrate/segmentation-fold
	 $ cd segmentation-fold
	 $ cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=~/.local .
	 $ make
	 $ make test
	 $ make install

You will then find the binary in the directory:
	<home directory>/.local/bin/

### USAGE ###

	Usage: ./segmentation-fold -h minimum_hairpin_length -s [SEQUENCE]
	Usage: ./segmentation-fold -h minimum_hairpin_length -f [FASTA_FILE]
	  * Note: If FASTA_FILE and SEQUENCE are not provided,
			  the program will read from STDIN.


	The following parameters can be used:
	  -h HAIRPINSIZE             Minimum hairpin size
	  -f FASTA_FILE              Location to FASTA_FILE
	  -s SEQUENCE                Specific RNA SEQUENCE (overrules -f)
	  -k [1/0]                   Enable/disable K-turn predictions
	  -m MOTIFS_FILE             Use custom "motifs.xml"-syntaxed file
	  -t [0/1]                   Use novel (=0) or original (=1) traceback
	  -v                         Shows the version and license


### BUGS ###

Please report encountelred problems at:

<https://github.com/yhoogstrate/segmentation-fold/issues>


### AUTHORS ###

All people that have contributed to the development of segmentation-fold are
mentioned in the '[AUTHORS](https://github.com/yhoogstrate/segmentation-fold/blob/master/AUTHORS)' file. You can read it by using on of the
following commands:

	 $ less AUTHORS

### DEVELOPER NOTES ###

We encourage users, researchers and programmers to contribute to this
free and open source project. This can be achieved by reporting bugs and
commiting code to this github repository. To streamline and archive
communication in an univocal way, we also encourge you to only use this
channel for contribution.

If any developer wants to contribute to segmentation-fold, please use
the following documentation standard:
 <http://www.stack.nl/~dimitri/doxygen/index.html>

All source-code is formatted using Asteric Style:
 <http://astyle.sourceforge.net/>
The corresponding configuration file is available under
'[src/.astylerc](https://github.com/yhoogstrate/segmentation-fold/tree/master/share/.astylerc)'.

For more developer notes you should consult [CONTRIBUTING.md](https://github.com/yhoogstrate/segmentation-fold/blob/master/CONTRIBUTING.md)
