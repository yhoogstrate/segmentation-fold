#!/usr/bin/env python

"""
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2016 Youri Hoogstrate

This file is part of segmentation-fold

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

import pysam,re,logging

class ExtractBoxedSequences:
    """
    """

    logger = logging.getLogger("segmentation_fold_utils::ExtractBoxedSequences")

    box1  = "box1"
    box2  = "box2"
    box1r = "box1r"
    box2r = "box2r"

    prog_actg_seq = re.compile("^[ACTGactg]+$")

    def __init__(self, fasta_input_file, bed_input_file, fasta_output_file, max_inner_dist, bp_extension):
        self.fasta_input_file = fasta_input_file
        self.bed_input_file = bed_input_file
        self.max_inner_dist = int(max_inner_dist)
        self.bp_extension = int(bp_extension)
    
    def reverse_complement(self,sequence):
        tt = {
            'A':'T',
            'a':'t',
            'C':'G',
            'c':'g',
            'T':'A',
            't':'a',
            'G':'C',
            'g':'c'
        }
        sequence_rc = ''
        for base in sequence:
            sequence_rc += tt[base]
        
        return sequence_rc[::-1]
    
    def insert_box1(self, lib, position):
        """
        Inserts box1 or box1r - this is usually the box that is lagging behind
        """
        lib.append(position)

    def update_using_box2(self, lib, position, fasta_output_file):
        """
        For the D-box, find whether there are C-boxes less than max-inner-dist away
        If so, link them

        If an iteration is beyond a certain cbox, get it out.
        If it has linked D-box, store it or print it
        """
        to_pop = []
        for i in range(len(lib)):
            item = lib[i]
            
            if item[0] != position[0] or item[2]+self.max_inner_dist < position[1]:
                if len(item[4]) != 0:
                    for box2 in item[4]:
                        _start = max(0,item[1]-self.bp_extension)
                        _end = box2[2]+self.bp_extension
                        seq = self.ref.fetch(item[0],_start,_end)
                        
                        if self.prog_actg_seq.match(seq):
                            if position[3] == "-":
                                seq = self.reverse_complement(seq)
                                fasta_output_file.write(">"+item[0]+":"+str(_end)+"-"+str(_start)+"(-)\n")
                            else:
                                fasta_output_file.write(">"+item[0]+":"+str(_start)+"-"+str(_end)+"(+)\n")
                            fasta_output_file.write(seq+"\n")
                        else:
                           self.logger.debug("Invalid sequencae was found: "+item[0]+":"+str(_start)+"-"+str(_end)+"\n"+seq)
                    to_pop.append(i)
            else:
                lib[i][4].append(position)
        
        to_pop.reverse()
        for i in to_pop:
            lib.pop(i)

    def run(self, fasta_output_file):
        self.ref = pysam.FastaFile(self.fasta_input_file)
        boxes_forward = []
        boxes_reverse = []
        
        self.bed_input_file.seek(0)
        for line in self.bed_input_file:
            line = line.strip()
            if len(line) > 0:
                params = line.split("\t")
                if len(params) > 1:
                    _chr, _start, _end, _name, _strand = params
                    _start = int(_start)
                    _end = int(_end)
                    _box, _seq = _name.split(":",1)
                    if _box == "box1-f":
                        boxes_forward.append([_chr,_start,_end,_strand,[]])
                    elif _box == "box2-f":
                        self.update_using_box2(boxes_forward,[_chr,_start,_end,_strand],fasta_output_file)
                    
                    elif _box == "box2-r":
                        boxes_reverse.append([_chr,_start,_end,_strand,[]])
                    elif _box == "box1-r":
                        self.update_using_box2(boxes_reverse,[_chr,_start,_end,_strand],fasta_output_file)
        
        self.update_using_box2(boxes_forward,["NOTEXIST",-1,-1,"+"],fasta_output_file)# Trigger last one to be exported as well
        self.update_using_box2(boxes_reverse,["NOTEXIST",-1,-1,"-"],fasta_output_file)# Trigger last one to be exported as well
