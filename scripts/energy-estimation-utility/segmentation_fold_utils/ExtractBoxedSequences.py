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

    prog_actg_seq = re.compile("^[ACTGUactgu]+$")

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

    def update_cboxes(self, lib, position, fasta_output_file):
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
                if item[4] != None:
                    _start = max(0,item[1]-self.bp_extension)
                    _end = item[4][2]+self.bp_extension
                    seq = self.ref.fetch(item[0],_start,_end)
                    if self.prog_actg_seq.match(seq):
                        fasta_output_file.write(">"+item[0]+":"+str(_start)+"-"+str(_end)+"\n")
                        fasta_output_file.write(seq+"\n")
                    else:
                       self.logger.debug("Invalid sequencae was found: "+item[0]+":"+str(_start)+"-"+str(_end)+"\n"+seq)
                to_pop.append(i)
            else:
                lib[i][4] = position
        
        to_pop.reverse()
        for i in to_pop:
            lib.pop(i)

    def run(self, fasta_output_file):
        self.ref = pysam.FastaFile(self.fasta_input_file)
        boxes_forward = []
        
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
                        print boxes_forward
                        self.insert_box1(boxes_forward,[_chr,_start,_end,_strand])
                        print " --- "
                        print boxes_forward
                        print "------------------------------------\n\n"
                    elif _box == "box2-f":
                        self.update_cboxes(boxes_forward,[_chr,_start,_end,_strand],fasta_output_file)
                    else:
                        print "["+_box+"]"
        
        self.update_cboxes(boxes_forward,["NOTEXIST",-1,-1,"+"],fasta_output_file)# Trigger last one to be exported as well




