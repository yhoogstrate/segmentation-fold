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



import random

from segmentation_fold_utils.RNA import RNA
from segmentation_fold_utils.FastaFile import FastaFile
from segmentation_fold_utils.BinarySplit import BinarySplit

from xml.dom.minidom import parseString


class XMLFile:
    def __init__(self,xml_input_file):
        self.xml_input_file = xml_input_file
        
        self.segments = {}
        #self.segmentloops = {}
        self.tests = []
        
        self.parse()
    
    def parse(self):
        data = self.xml_input_file.read()
        
        root = parseString(data)
        segments = root.getElementsByTagName('segments')[0]
        
        for segment in segments.getElementsByTagName('segment'):        # Read segments
            name        = segment.getElementsByTagName('id')[0].firstChild.data
            p5_sequence = segment.getElementsByTagName('sequence_5prime')[0].firstChild.data
            p3_sequence = segment.getElementsByTagName('sequence_3prime')[0].firstChild.data
            bonds       = segment.getElementsByTagName('bonds')[0].firstChild.data
            energy      = segment.getElementsByTagName('energy')[0].firstChild.data
            
            self.segments[name] = {'5':p5_sequence,'3':p3_sequence,'x':bonds,'e':energy}
        
        for rna_xml in root.getElementsByTagName('rna'):                # Read RNAs
            name            = rna_xml.getElementsByTagName('title')[0].firstChild.data
            sequence        = rna_xml.getElementsByTagName('sequence')[0].firstChild.data
            
            organism_xml    = rna_xml.getElementsByTagName('organism')
            structures_xml  = rna_xml.getElementsByTagName('structures')[0]
            
            if(organism_xml and organism_xml[0].firstChild):
                organism = organism_xml[0].firstChild.data
            else:
                organism = None
            
            structures = []
            for structure_xml in structures_xml.getElementsByTagName('structure'):
                structures.append(self.parse_structure(structure_xml))
            
            self.tests.append(RNA(name,sequence,organism,structures))
        
        return True
    
    def parse_structure(self,structure_xml):
        structure = {'associated_segments':[],'dot_bracket':None}
        
        associated_segments = structure_xml.getElementsByTagName('associated_segments')[0]
        for associated_segment in associated_segments.getElementsByTagName('link'):
            structure['associated_segments'].append(associated_segment.firstChild.data)
        
        dot_bracket = structure_xml.getElementsByTagName('dot_bracket')
        if(dot_bracket and dot_bracket[0].firstChild):
            structure['dot_bracket'] = dot_bracket[-1].firstChild.data
        else:
            structure['dot_bracket'] = None
        
        return structure
    
    def __iter__(self):
        for i in range(len(self.tests)):
            yield self.tests[i]
    
    def sanitize_sequence(self,seq,rev=False):
        seq = seq.strip().upper()
        seq = seq.replace("T","U")
        seq = seq.replace(" ","")
        if rev:
            seq = seq[::-1]
        
        return seq
    
    def sanitize_name(self,name):
        name = name.replace(" ","_").replace("\t","_")
        
        charset = '-_.abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
        return ''.join(c for c in name if c in charset)
        
    def shuffle_sequence(self,sequence,segments):
        """Should preserve all subsequences of the motif
        """
        # Find all subsequences that are part of a segment
        subsequences = []
        for segment in segments.values():
            subsequences.append(self.sanitize_sequence(segment['5']))
            subsequences.append(self.sanitize_sequence(segment['3'],True))
        
        # Order those subsequence on size, so that searching will first find the largest subsequence in the sequence
        subsequences.sort(key=len,reverse=True)
        
        # For every position in the sequence, check whether this equals one of the subsequences ordered on size
        slices = []
        i = 0
        while i < len(sequence):
            add = True
            j = 0
            while(add and j < len(subsequences)):
                subsequence = subsequences[j]
                if sequence[i:i+len(subsequence)] == subsequence:
                    add = False
                
                j += 1
            
            # If it not found, add the base
            if add:
                slices.append(sequence[i].lower())
                i += 1
            else:# Otherwise, add the subsequence as block
                slices.append(subsequence)
                i += len(subsequence)
        
        random.shuffle(slices)
        return "".join(slices)

    def get_combinations(self,fasta_input_file):
        if fasta_input_file == None:
            for item in self:
                for segment in item.get_unique_associated_segments():
                    yield item.name, item.sequence, segment, self.segments[segment]
        else:
            fasta = FastaFile(fasta_input_file)
            for sequence in fasta:
                for segment in self.segments.keys():
                    yield sequence['name'], sequence['sequence'], segment, self.segments[segment]
    
    def estimate_energy(self,temp_dir,segmentation_fold,threads,precision,randomize,sequences_from_fasta_file,dbn_output_file):
        for sequence_name,sequence,segment_name,segment in self.get_combinations(sequences_from_fasta_file):
            sequence_name_s = self.sanitize_name(sequence_name)
            segment_name_s = self.sanitize_name(segment_name)
            
            if randomize <= 0:
                dbn_output_file.write(">"+sequence_name+" x "+segment_name+"\n")
                dbn_output_file.write(sequence+"\n")
                
                b = BinarySplit(
                    segmentation_fold,
                    temp_dir+"/segments_"+sequence_name_s+"_"+segment_name_s+"_",
                    sequence,
                    {segment_name:segment},
                    precision, (3.5/2), threads) 
                for transition in b.find_transitions():
                    dbn_output_file.write(transition['structure_max']+"\t"+transition['structure_min']+"\t"+str(transition['energy'])+"\n")
            else:
                for i in range(1,randomize+1):
                    sequence_r = self.shuffle_sequence(sequence,{segment_name:segment})
                    sequence_name_r = self.sanitize_name(sequence_name+" (shuffle iteration: "+str(i)+")")
                    
                    dbn_output_file.write(">"+sequence_name_r+") x "+segment_name+"\n")
                    dbn_output_file.write(sequence_r+"\n")
                    
                    b = BinarySplit(
                        segmentation_fold,
                        temp_dir+"/segments_"+self.sanitize_name(sequence_name_r)+"_"+segment_name_s+"_",
                        sequence_r,
                        {segment_name:segment},
                        precision, (3.5/2), threads) 
                    for transition in b.find_transitions():
                        dbn_output_file.write(transition['structure_max']+"\t"+transition['structure_min']+"\t"+str(transition['energy'])+"\n")
