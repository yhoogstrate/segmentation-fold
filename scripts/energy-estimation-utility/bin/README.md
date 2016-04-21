energy-estimation-utilityby 
-------------------------
This tool will try to figurea out for a givin segment within a given sequence, the (minimal/maximal) energy values that is necessary to become part of the optimal overall 2D structure. This limit is not by definition the energy value that corresponds to the actual segment, but should be - on the assumption that the other parameters are correct - but sequences in which the energy is ... should not contain this substructre whilst containing the subsequence.

As example we have the following two miRNAs:
 - * from analysis
```
5') AC-d-box-GACUGUACCUGACA
    \\       |||| || ||||  C
3')  UG-c-boxCUGAGAUAGACUAA
```


scan-for-segments
-----------------
This tool allows to search for the presence of segments or segmentloops within any arbitrary fasta sequence. It can be set in the *in-depth* mode which more or less only tries to find the possibilitiy, for each given segment, whether it can be folded if the total energy contribution of the segment would be minus inifinty. In the *in-depth* method it will 

@ todo change to 'in-general' and 'energy-restricted'
