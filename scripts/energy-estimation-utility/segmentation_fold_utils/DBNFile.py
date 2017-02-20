#!/usr/bin/env python

"""
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2016 Youri Hoogstrate

This file is part of segmentation-fold and originally taken from
yh-kt-fold.

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


import random,re,logging,subprocess
import pysam

from segmentation_fold_utils.GeneAnnotation import GeneAnnotation


class DBNFile:
    
    logger = logging.getLogger("segmentation-fold-utils::DBNFile")
    
    def __init__(self,dbn_input_file,estimate_energy_results):
        """
        if 'estimate_energy_results' is false, the (multi-)dbn format
        should be:

        The Dot-Bracket Notation format (.dbn) is a commonly used format to describe 2D RNA structures, in particular as output for structure predictors that can't predict pseudoknots. Every entry consists of three lines:

        1. The header line, to provide information about the sequence (e.g. name, free energy).
        2. The sequence line, to provide the sequence of which the structures will be described.
        3. The structure (Dot-Bracket) line, describing the structure of the sequence.
        Two nucleotides that form a bond are indicated with encloding parenthesis and dots represent unbound nucleotides.
        Pseudoknots are indicated with { .. } or [ .. ].
        
        else if 'estimate_energy_results' is true, the format should be:
         - first line should be a fasta header line, starting with '>'
         - second line should be the corresponding sequence
         - the following, optional, lines should be:
            * Structure-string-min\tStructure-string-max\tdelta-G
        
        It could thus very well be that there is no structure string for
        a particular sequence because there were no differerences in the
        overal structure detected by changing the energy value of the
        segment
        """
        
        self.dbn_input_file = dbn_input_file
        self.estimate_energy_results = estimate_energy_results
    
    def parse(self):
        self.dbn_input_file.seek(0)
        n_item  = None
        
        for line in self.dbn_input_file:
            line = line.strip()
            if len(line) > 0:
                if line[0] == '>':
                    if n_item != None:
                        yield n_item
                    
                    if self.estimate_energy_results:
                        n_item = {'name':line,'sequence':None,'transitions':[]}
                    else:
                        n_item = {'name':line,'sequence':None,'sequence':None}
                else:
                    if n_item['sequence'] == None:
                        n_item['sequence'] = line
                    else:
                        if self.estimate_energy_results:
                            params = line.split("\t")
                            params[2] = float(params[2])
                            n_item['transitions'].append(params)
                        else:
                            n_item['sequence'] = line
        
        if n_item != None:
            yield n_item
    
    def annotate_read_counts(self,bam_input_file,regex,dbn_output_file):
        regex_prog = re.compile(regex)
        self.check_bam_index(bam_input_file)
        alignment = pysam.Samfile(bam_input_file)
        
        k = 0
        l = 0
        for structure in self.parse():
            k += 1
            j = 0
            
            for position_g in regex_prog.finditer(structure['name']):
                position = [position_g.group(1),int(position_g.group(2)),int(position_g.group(3))]
                if position[1] > position[2]:
                    position = [position[0], position[2], position[1]]
                
                for read in alignment.fetch(position[0], position[1], position[2]):
                    j += 1
                
            dbn_output_file.write(structure['name']+' (aligned reads '+bam_input_file+': '+str(j)+')'+"\n")
            dbn_output_file.write(structure['sequence']+"\n")
            for transition in structure['transitions']:
                transition_f = transition
                transition_f[2] = str(transition_f[2])
                dbn_output_file.write("\t".join(transition_f)+"\n")
                l += j
        
        if k == 0:# Click does not produce a file if nothing is written to it, not even with open(...,'w+)
            dbn_output_file.write("")
        
        self.logger.info("Written "+str(k)+" sequences with "+str(l)+" annotated reads to "+dbn_output_file.name)
    
    def check_bam_index(self,bam_file):
        fh = pysam.Samfile(bam_file)
        
        if(len(fh.references) > 0):
            try:
                fh.fetch(fh.references[0], 0, 0)
                fh.close()
            except:
                fh.close()
                self.logger.warning('Missing index for bam file, trying to index with samtools')
                try:
                    pysam.index(bam_file)
                except:
                    self.logger.warning('Couldn\'t indexing BAM file with samtools: '+bam_file+'. Are you sure samtools is installed?')
        else:
            self.logger.warning('Missing chromosomes in BAM file - is it empty?')

    def filter_annotated_entries(self,regex,bed_input_file,dbn_output_file_o,dbn_output_file_n):
        genes = GeneAnnotation()
        genes.parse_bed(bed_input_file)
        
        i = 0
        j = 0
        
        regex_prog = re.compile(regex)
        
        for structure in self.parse():
            annotations = set([])
            for position_g in regex_prog.finditer(structure['name']):
                position = [position_g.group(1),int(position_g.group(2)),int(position_g.group(3))]
                if position[1] > position[2]:
                    position = [position[0], position[2], position[1]]
                
                annotations = annotations.union(genes.get_annotations(position[0],position[1],position[2]))
            
            annotations = list(annotations)
            if len(annotations) == 0:
                i += 1
                dbn_output_file_n.write(structure['name']+"\n")
                dbn_output_file_n.write(structure['sequence']+"\n")
                for transition in structure['transitions']:
                    transition_f = transition
                    transition_f[2] = str(transition_f[2])
                    dbn_output_file_n.write("\t".join(transition_f)+"\n")

            else:
                j += 1
                dbn_output_file_o.write(structure['name']+' (overlap in '+bed_input_file.name+': '+",".join(annotations)+')'+"\n")
                dbn_output_file_o.write(structure['sequence']+"\n")
                for transition in structure['transitions']:
                    transition_f = transition
                    transition_f[2] = str(transition_f[2])
                    dbn_output_file_o.write("\t".join(transition_f)+"\n")
        
        # Click does not produce a file if nothing is written to it, not even with open(...,'w+)
        if i == 0:
            dbn_output_file_n.write("")
        if j == 0:
            dbn_output_file_o.write("")
        
        self.logger.info("Written "+str(i)+" (unannotated) entries to "+dbn_output_file_n.name+" and "+str(j)+" (annotated) entries to "+dbn_output_file_o.name)

    def filter_by_energy(self,dbn_output_file_larger_or_equal,dbn_output_file_smaller,energy):
        l = 0
        s = 0
        
        for structure in self.parse():
            tl = 0
            ts = 0
            for transition in structure['transitions']:
                transition_s = "\t".join([str(x) for x in transition])
                
                if transition[2] >= energy:
                    if tl == 0:
                        # print headers
                        dbn_output_file_larger_or_equal.write(structure['name']+"\n")
                        dbn_output_file_larger_or_equal.write(structure['sequence']+"\n")
                    
                    dbn_output_file_larger_or_equal.write(transition_s+"\n")
                    tl += 1
                else:
                    if ts == 0:
                        # print headers
                        dbn_output_file_smaller.write(structure['name']+"\n")
                        dbn_output_file_smaller.write(structure['sequence']+"\n")
                    
                    dbn_output_file_smaller.write(transition_s+"\n")
                    ts += 1
                
            if ts+tl == 0:
                dbn_output_file_smaller.write(structure['name']+"\n")
                dbn_output_file_smaller.write(structure['sequence']+"\n")
                pass
            
            l += tl
            s += ts
        
        # Click does not produce a file if nothing is written to it, not even with open(...,'w+)
        if l == 0:
            dbn_output_file_larger_or_equal.write("")
        if s == 0:
            dbn_output_file_smaller.write("")
        
        self.logger.info("Written "+str(l)+" (larger or equal than) entries to "+dbn_output_file_larger_or_equal.name+" and "+str(s)+" (smaller than) entries to "+dbn_output_file_smaller.name)
