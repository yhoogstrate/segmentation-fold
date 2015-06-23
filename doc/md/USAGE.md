
# USAGE #

Usage: segmentation-fold -s [SEQUENCE]
Usage: segmentation-fold -f [FASTA_FILE]
   * Note: If FASTA_FILE and SEQUENCE are not provided,
           the program will read from STDIN.


The following parameters can be used:
  -s SEQUENCE                Specific RNA SEQUENCE (overrules -f)
  -f FASTA_FILE              Location to FASTA_FILE that contains the sequence

  -p                  [1/0]  Enable/disable segment prediction functionality

  -h HAIRPINSIZE        [3]  Minimum hairpin size, 0 or larger, default 3
  -x SEGMENTS_XML_FILE       Use custom  "segments.xml"-syntaxed file

  -V --version               Shows the version and license


If you encounter problems with this software, please send bug-reports to:
   <https://github.com/yhoogstrate/segmentation-fold/issues>

## galaxy ##

segmentation-fold will become available as galaxy tool.
