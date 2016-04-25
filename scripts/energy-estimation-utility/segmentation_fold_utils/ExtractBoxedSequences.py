#!/usr/bin/env python

"""
license: GPL3
"""

import pysam,re

class ExtractBoxedSequences:
    """
    """

    logger = logging.getLogger("segmentation_fold_utils::ExtractBoxedSequences")

    box1  = "box1"
    box2  = "box2"
    box1r = "box1r"
    box2r = "box2r"
    
    def __init__(self, fasta_input_file, bed_input_file, fasta_output_file, max_inner_dist, bp_extension):
        self.inner_dist = inner_dist

    def insert_box1(self, lib, position):
    """
    Inserts box1 or box1r - this is usually the box that is lagging behind
    """
        position.append(None)
        previous_cboxes.append(position)

    def update_cboxes(lib, position):
        """
        For the D-box, find whether there are C-boxes less than max-inner-dist away
        If so, link them

        If an iteration is beyond a certain cbox, get it out.
        If it has linked D-box, store it or print it
        """
        to_pop = []
        for i in range(len(lib)):
                item = lib[i]
                
                if item[0] != position[0] or item[1]+len(motif_a)+max_inner_dist < position[1]:
                        #print "Cleaning",item
                        if item[2] != None:
                                #print "C/D box at ",item," d=",item[2][1]-item[1]-len(motif_a)
                                seq = ref.fetch(item[0],item[1]-bp_extension,item[2][1]+bp_extension+len(motif_b)).replace("T","U").replace("t","u")
                                if re.match("^[ACTGUactgu]+$",seq):
                                        print ">"+item[0]+":"+str(item[1]-bp_extension)+"-"+str(item[2][1]+bp_extension+len(motif_b)-1)
                                        print seq
                        
                        to_pop.append(i)
                else:
                        lib[i][2] = position
        
        to_pop.reverse()
        for i in to_pop:
                lib.pop(i)


    def run(self, output_fasta_file):
        ref = pysam.FastaFile(self.fasta_input_file)

            # Seems obsolete:
            #motif_a = 'RUGAUG'.replace('U','T')
            #motif_b = 'CUGA'.replace('U','T')
            #bp_extension = 10

            previous_cboxes = []

                with open(file_fwd,"r") as fh:
                        for line in fh:
                                line = line.strip()
                                if len(line) > 0:
                                        params = line.split("\t")
                                        if len(params) > 1:
                                                params[0] = params[0].split(":")
                                                params[0][1] = int(params[0][1])
                                                
                                                if params[1] == "C-box":
                                                        insert_box1(previous_cboxes,params[0])
                                                
                                                elif params[1] == "D-box":
                                                        update_cboxes(previous_cboxes,params[0])




