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


class RNA:
    def __init__(self,name,sequence,organism,structures):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.structures = structures
    
    def get_sequence(self):
        return self.sequence
    
    def get_structures(self):
        return self.structures
    
    def get_unique_associated_segments(self):
        segments = []
        
        for structure in self.structures:
            for associated_segment in structure['associated_segments']:
                segments.append(associated_segment)
        
        return list(set(segments))
        
