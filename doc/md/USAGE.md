
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
