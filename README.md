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

The tool segmentation-fold is a program that does RNA 2D-structure
prediction using a modified version of the Zuker algorithm. This
modification includes a new "structure element" that is capable of
folding a certain type of motifs (predefined sequences).
	The program includes some motifs; the K-turns. The thermodynamic
parameters (free Gibbs energy levels) have been estimated using a
biassed computational approach and lack accuracy.


### INSTALLATION ###

Please follow the instructions given in the '[INSTALL](https://github.com/yhoogstrate/segmentation-fold/blob/master/INSTALL)' file, e.g. by
using one of the following commands:

	 $ vi INSTALL
	 $ less INSTALL
	 $ cat INSTALL
	 $ pico INSTALL
	 $ nano INSTALL


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

Please report encountered problems at:

<https://github.com/yhoogstrate/segmentation-fold/issues>


### AUTHORS ###

All people that have contributed to the development of segmentation-fold are
mentioned in the '[AUTHORS](https://github.com/yhoogstrate/segmentation-fold/blob/master/AUTHORS)' file. You can read it by using on of the
following commands:

	 $ vi AUTHORS
	 $ less AUTHORS
	 $ cat AUTHORS
	 $ pico AUTHORS
	 $ nano AUTHORS


### DEVELOPER NOTES ###

If any developer wants to contribute to segmentation-fold, please use
the following documentation standard:
 <http://www.stack.nl/~dimitri/doxygen/index.html>

All source-code is formatted using Asteric Style:
 <http://astyle.sourceforge.net/>
The corresponding configuration file is available under
'[src/.astylerc](https://github.com/yhoogstrate/segmentation-fold/tree/master/share/.astylerc)'.
