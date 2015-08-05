#!/usr/bin/env python

"""
@file scripts/energy-estimation-utility/energy_estimation_utility/DataController.py

@date 2015-07-20

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



from xml.dom.minidom import parseString
from energy_estimation_utility.RNA import *



class DataController:
    def __init__(self,xmlfile):
        self.xmlfile = xmlfile
        
        self.segments = {}
        self.tests = []
        
        self.parse()
    
    def parse(self):
        fh = open(self.xmlfile,'r')
        data = fh.read()
        fh.close()
        
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