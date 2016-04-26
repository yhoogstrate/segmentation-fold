#!/usr/bin/env python

"""
license: GPL3
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

    def __init__(self, fasta_input_file, bed_input_file, fasta_output_file, max_inner_dist, bp_extension):
        self.fasta_input_file = fasta_input_file
        self.bed_input_file = bed_input_file
        self.max_inner_dist = int(max_inner_dist)
        self.bp_extension = int(bp_extension)

    def insert_box1(self, lib, position):
        """
        Inserts box1 or box1r - this is usually the box that is lagging behind
        """
        position.append(None)
        lib.append(position)

    def update_cboxes(self, lib, position):
        """
        For the D-box, find whether there are C-boxes less than max-inner-dist away
        If so, link them

        If an iteration is beyond a certain cbox, get it out.
        If it has linked D-box, store it or print it
        """
        to_pop = []
        for i in range(len(lib)):
            item = lib[i]
            print item,"=box=",position
            
            if item[0] != position[0] or item[2]+self.max_inner_dist < position[1]:
                if item[4] != None:
                    print "**",item
                    #seq = self.ref.fetch(item[0],item[1]-self.bp_extension,position[2]+self.bp_extension)
                    #print seq
                    #if re.match("^[ACTGUactgu]+$",seq):
                    #    print ">"+item[0]+":"+str(item[1]-self.bp_extension)+"-"+str(item[2][1]+self.bp_extension+len(motif_b)-1)
                    #    print seq
                    #else:
                     #   self.logger.debug("Invalid sequencae was found: "+position[0]+":"+str(position[1])+"-"+str(position[2])+"\n"+seq)
                to_pop.append(i)
            else:
                lib[i][4] = position
        
        to_pop.reverse()
        for i in to_pop:
            lib.pop(i)

    def run(self, output_fasta_file):
        self.ref = pysam.FastaFile(self.fasta_input_file)
        boxes_foward = []
        
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
                        self.insert_box1(boxes_foward,[_chr,_start,_end,_strand])
                    elif _box == "box1-r":
                        self.update_cboxes(boxes_foward,[_chr,_start,_end,_strand])
                    




