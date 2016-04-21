#!/usr/bin/env python

"""
@section LICENSE
<PRE>
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
</PRE>
"""


import logging
import pysam


class FindCDboxes:
	"""detecting C boxes in chr22 takes ~7min
	"""
	
	logger = logging.getLogger("segmentation_fold_utils::FindCDboxes")
	
	def __init__(self,genome,box1,box2,search_fwd,search_rev,inner_dist):
		self.genome = genome
		self.box1 = box1.upper().replace('U','T').strip('N')
		self.box2 = box2.upper().replace('U','T').strip('N')
		self.k = len(self.box1)
		self.l = len(self.box2)
		self.m = max(self.k, self.l)
		self.search_fwd = search_fwd
		self.search_rev = search_rev
		self.inner_dist = inner_dist
		
	
	def match_base(self,base_pattern, base_query):
		if base_query == "A":
			return base_pattern in "ARMWDHV"#+N
		elif base_query == "G":
			return base_pattern in "GRKSBDV"#+N
		elif base_query == "C":
			return base_pattern in "CYMSBHV"#+N
		elif base_query == "T":
			return base_pattern in "TYKWBDH"#+N
		else:
			return False
	
	def match(self,pattern, query):
		for i in range(len(pattern)):
			if not self.match_base(pattern[i], query[i]):
				return False
		return True
	
	def run(self,output_file):
		fh = open(output_file,"w")
		ref = pysam.FastaFile(self.genome)
		for chromosome in ref.references:
			# Look fwd
			n = ref.get_reference_length(chromosome)
			self.logger.debug(chromosome+": ("+str(n)+")")
			
			#k=7
			#n=10
			#i[0] = 0,7
			#i[0] = 1,8
			#i[0] = 2,9
			#i[0] = 3,10
			
			#n = 10
			
			for i in range(n-self.m+1):
				chunk1 = ref.fetch(chromosome,i,i+self.k)
				if self.match(self.box1,chunk1):
					fh.write(chromosome+":"+str(i)+"\tbox1\n")
				
				chunk2 = ref.fetch(chromosome,i,i+self.l)
				if self.match(self.box2,chunk2):
					fh.write(chromosome+":"+str(i)+"\tbox2\n")
			
			#print
			#print match(motif_a,"A"+"TGATG"),"- should be True"
			#print match(motif_a,"G"+"TGATG"),"- should be True"
			#print match(motif_a,"C"+"TGATG"),"- should be False"
			#print match(motif_a,"T"+"TGATG"),"- should be False"
		fh.close()
