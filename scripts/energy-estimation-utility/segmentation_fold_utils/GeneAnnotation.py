#!/usr/bin/env python

"""
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2016 Youri Hoogstrate

This file is part of segmentation-fold.

segmentation-fold is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

segmentation-fold is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""


import HTSeq
import logging

class GeneAnnotation:
    """Gene annotation is a virtual reference genome. It's only being
    used the map Genes to in order to index them quickly. 
    """
    logger = logging.getLogger("segmentation-fold-utils::GeneAnnotation")
    
    def __init__(self):
        self.n = 0
        self.gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    
    def add_annotation(self,gene,chromosome,start,stop):
        self.gas[HTSeq.GenomicInterval(chromosome,start,stop)] += gene
        self.n += 1
    
    def get_annotations(self,chromosome,start,end):
        entries = set()
        for annotation in self.gas[HTSeq.GenomicInterval(chromosome,start,end,'.')]:
            entries = entries.union(annotation)
        
        return entries
    
    def __len__(self):
        return self.n
    
    def __iter__(self):
        for chromosome_name,chromosome_obj in self.gas.chrom_vectors.items():
            for gene in list(reduce(lambda s1, s2: s1 | s2, [x[1] for x in self.gas[HTSeq.GenomicInterval(chromosome_name,0,chromosome_obj['.'].iv.end)].steps()])):
                yield gene
    
    def parse_bed(self,bed_input_file):
        bed_input_file.seek(0)
        
        for line in bed_input_file:
            line = line.strip()
            if len(line) > 0 and line[0] != "#":
                params = line.split("\t")
                if len(params) >= 6:
                    _chr = params[0]
                    _start = int(params[1])
                    _end = int(params[2])
                    _name = params[3]
                    _score = int(params[4])
                    _strand = params[5]
                else:
                    self.logger.warning("Invalid BED file - not enough columns")
                
                self.add_annotation(_name, _chr, _start, _end)
