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
two new "structure elements", the *segment* and *segmentloop* which are
capable of folding pre-defined substructure containing **multiple**
canonical and/or **non-canonical** pairings.

This allows predicting structures containing sub-structures like K-turns
and loop-E-motifs, which are also part of the implemented free energy
tables [segments.xml](https://github.com/yhoogstrate/segmentation-fold/blob/master/share/segmentation-fold/segments.xml).
These parameters (free Gibbs energy) have been estimated using in silico
predictions but lack accuracy due to the limited number of available
datapoints. In order to add novel structures to segmentation-fold, the
*segments.xml* file needs to be modified.

We have made a set of utilities using segmentation-fold available at
[./scripts/energy-estimation-utilities/](https://github.com/yhoogstrate/segmentation-fold/blob/master/scripts/energy-estimation-utility/).
Details on the utilities can be found at [https://github.com/yhoogstrate/segmentation-fold/blob/master/scripts/energy-estimation-utility/README.md](https://github.com/yhoogstrate/segmentation-fold/blob/master/scripts/energy-estimation-utility/README.md).
