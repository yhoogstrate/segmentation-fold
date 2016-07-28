
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

segmentation-fold with utilities is available for galaxy at the following url:
[https://toolshed.g2.bx.psu.edu/view/yhoogstrate/segmentation_fold/](https://toolshed.g2.bx.psu.edu/view/yhoogstrate/segmentation_fold/),
and the code is maintained at:
[https://github.com/ErasmusMC-Bioinformatics/segmentation_fold_galaxy_wrapper](https://github.com/ErasmusMC-Bioinformatics/segmentation_fold_galaxy_wrapper)
