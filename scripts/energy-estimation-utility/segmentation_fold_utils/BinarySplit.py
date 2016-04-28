#!/usr/bin/env python

"""
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2016 Youri Hoogstrate

This file is part of segmentation-fold.

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



from energy_estimation_utility import FoldController
from energy_estimation_utility.FoldController import *



class BinarySplit:
    def __init__(self,binary,xml_file_prefix,sequence,arg_associated_segments,precision=0.005,max_per_base=(3.5/2),threads=1):
        self.binary = binary
        self.xml_file_prefix = xml_file_prefix
        
        self.sequence = sequence
        self.associated_segments = arg_associated_segments
        
        self.precision = precision
        self.max_per_base = max_per_base
        
        self.min_energy = - abs(max_per_base) * len(sequence)
        
        self.threads_per_instance = threads
    
    def update_energy_in_segments(self,energy):
        segments = {}
        for name, data in self.associated_segments.items():
            segments[name] = self.associated_segments[name]
            segments[name]['e'] = energy
        
        return segments
    
    def find_transitions(self,min_energy = {'energy':None,'results':None},max_energy = {'energy':None,'results':None},stats = {'state':None,'recursion_depth':1}):
        if(min_energy['energy'] == None):
            min_energy['energy'] = self.min_energy
            max_energy['energy'] = abs(self.min_energy)
        
        precision = max_energy['energy'] - min_energy['energy']
        
        xml_file = self.xml_file_prefix + str(stats['recursion_depth']) + "_" + (stats['state'] if stats['state'] != None else "min")+".xml"
        min_associated_segments = self.update_energy_in_segments(min_energy['energy'])
        fc_min = FoldController(self.binary,xml_file,self.sequence,min_associated_segments,self.threads_per_instance)
        min_energy['results'] = fc_min.fold()
    
        xml_file = self.xml_file_prefix + str(stats['recursion_depth']) + "_" + (stats['state'] if stats['state'] != None else "max")+".xml"
        max_associated_segments = self.update_energy_in_segments(max_energy['energy'])
        fc_max = FoldController(self.binary,xml_file,self.sequence,max_associated_segments,self.threads_per_instance)
        max_energy['results'] = fc_max.fold()
    
        #if min_energy['results']['dot_bracket'] != max_energy['results']['dot_bracket']:
        if min_energy['results']['number_segments'] != max_energy['results']['number_segments']:
            splitpoint = (max_energy['energy'] + min_energy['energy']) / 2
            
            if(precision >= self.precision):
                results_min = self.find_transitions(min_energy, {'energy':splitpoint,'results':None},{'state':'min','recursion_depth':stats['recursion_depth']+1})
                results_max = self.find_transitions({'energy':splitpoint,'results':None}, max_energy,{'state':'max','recursion_depth':stats['recursion_depth']+1})
                
                return results_min + results_max
            else:
                
                return [{'structure_min':min_energy['results']['dot_bracket'],'structure_max':max_energy['results']['dot_bracket'],'energy':splitpoint}]
        else:
            return []
