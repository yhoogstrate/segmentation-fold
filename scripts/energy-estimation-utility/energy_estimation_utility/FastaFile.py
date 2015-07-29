#!/usr/bin/env python

"""
@file scripts/energy-estimation/energy_estimation_utility/FastaFile.py

@date 2015-07-29

@author Youri Hoogstrate

@section LICENSE
<PRE>
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2015 Youri Hoogstrate

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
</PRE>
"""

class FastaFile:
    def __init__(self,filename):
        self.filename = filename
    
    def parse(self):
        sequence = None
        
        with open(self.filename,'r') as self.fh:
            for line in self.fh:
                line_s = line.strip()
                if(len(line_s) > 0):
                    if(line[0] == '>'):
                        if(sequence):
                            yield sequence
                        sequence = {'name':line_s[1:],'sequence':''}
                    else:
                        sequence['sequence'] += line_s.upper()
        
        if(sequence):
            yield sequence
    
    def __iter__(self):
        for sequence in self.parse():
            yield sequence
