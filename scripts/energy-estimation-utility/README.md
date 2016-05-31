segmentation-fold-utils
=======================
The computer program `segmentation-fold-utils` contains the following sub programs that allow downstream analysis using segmentation-fold:

1. **find-boxes**                Finds all occurances of two given boxes (sequence motifs) within a FASTA file
2. **extract-boxed-sequences**   Extracts boxes sequences from `bed_input_file` which has to be created with 'find-box', part of this utility
3. **estimate-energy**           Estimates whether a certain Segment(Loop) is present and for which delta-G this transistion takes place
4. **add-read-counts**           Annotate sequences by adding the read counts from a bam file, within a region contained in the fasta header of the dbn file
5. **filter-annotated-entries**  Split entries into two files based on whether they overlap annotations in a bed file

find-boxes
----------
*Finds all occurances of two given boxes (sequence motifs) within a FASTA file*

This tool allows the user to scan through the entire genome and enlist all
occurances of a certain box (default C- and D-box) and reports them into
a BED file.

**input**

The provided boxes must be in concorance with the FASTA format:
<https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation>

**output**

The results are in BED format, meaning: "0-based, half open". Hence, the
first base in the genome is detanoted as: `chr1 \t 0 \t 1`.

extact-boxed-sequences
----------------------
*Extracts boxes sequences from `bed_input_file` which has to be created with 'find-box', also part of this utility*

The user can use this utility to extract sequences containing the boxes provided in the bed file by `find-boxes`.

**input**

Important information about the input:

 - `FASTA_INPUT_FILE` can be any generic FASTA file that can be read with pysam. This means that if the sequence is split into multiple lines, they must all be at the same length.
 - `BED_INPUT_FILE` the bed file should be provided by `find-boxes` as it properly denotes the names (box1-f, box1-r, box2-f and box2-r) which are used for extraction.
 - `-d, --max-inner-dist INTEGER` Only sequences for which the distance in bases between the boxes is smaller than this distance, will be extracted. Boxes are excluded from this distance.
 - `-e, --bp-extension INTEGER` Each sequence will be exteded with:
  * The boxes
  * An optional number of bases provided with this argument

**output**

Be aware that there can be overlapping sequences. For example, if you started box1=`TTTT` and box2=`CCCC` with the following sequence, you will extract 2 sequences:

```>seq
gagagaTTTTgagagaTTTTgagagagagagagagaCCCCgaga
```

Namely:

```TTTTgagagaTTTTgagagagagagagagaCCCC
```

and

```          TTTTgagagagagagagagaCCCC
```

estimate-energy
---------------
*Estimates whether a certain Segment(Loop) is present and for which delta-G this transistion takes place*

This tool will try to figurea out for a givin segment within a given sequence, the (minimal/maximal) energy values that is necessary to become part of the optimal overall 2D structure. This limit is not by definition the energy value that corresponds to the actual segment, but should be - on the assumption that the other parameters are correct - but sequences in which the energy is ... should not contain this substructre whilst containing the subsequence.

As example we have the following two miRNAs:
 - * from analysis
```
5') AC-d-box-GACUGUACCUGACA
    \\       |||| || ||||  C
3')  UG-c-boxCUGAGAUAGACUAA
```


This tool allows to search for the presence of segments or segmentloops within any arbitrary fasta sequence. It can be set in the *in-depth* mode which more or less only tries to find the possibilitiy, for each given segment, whether it can be folded if the total energy contribution of the segment would be minus inifinty. In the *in-depth* method it will 

@ todo change to 'in-general' and 'energy-restricted'


add-read-counts
---------------
*Annotate sequences by adding the read counts from a bam file, within a region contained in the fasta header of the dbn file*



filter-annotated-entries
------------------------
*Split entries into two files based on whether they overlap annotations in a bed file*


