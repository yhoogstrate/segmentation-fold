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

import re,logging

class FastaFile:
    logger = logging.getLogger("segmentation-fold-utils::FastaFile")
    prog = re.compile("^>.*? [^ ]")
    
    def __init__(self,filename):
        self.filename = filename
    
    def fix_fasta_file(self, new_file_handle):
        diffs = 0
        with open(self.filename,'r') as fh:
            for line in fh:
                if self.prog.match(line):
                    new_file_handle.write(line.replace(" ","_"))
                    diffs += 1
                else:
                    new_file_handle.write(line)
        
        return diffs
    
    def parse(self):
        sequence = None
        
        with open(self.filename,'r') as fh:
            for line in fh:
                line_s = line.strip()
                if(len(line_s) > 0):
                    if(line[0] == '>'):
                        if(sequence):
                            if re.match("^[ACTGUactug]+$",sequence['sequence']):
                                yield sequence
                        sequence = {'name':line_s[1:],'sequence':''}
                    else:
                        sequence_line = line_s.upper().replace("\n","").strip()
                        sequence['sequence'] += sequence_line
        
        if(sequence):
            if re.match("^[ACTGUactug]+$",sequence['sequence']):
                yield sequence
    
    def __iter__(self):
        for sequence in self.parse():
            yield sequence
