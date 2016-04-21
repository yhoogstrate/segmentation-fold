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


import pysam


class FindCDboxes:
	"""detecting C boxes in chr22 takes ~7min
	"""
	def __init__(self,genome,box1,box2,search_fwd,search_rev,inner_dist):
		pass
	
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
			if not match_base(pattern[i], query[i]):
				return False
		return True
	
	def run(self):
		ref = pysam.FastaFile(fa_file)
		for chromosome in ref.references:
			# Look fwd
			n = ref.get_reference_length(chromosome)
			print chromosome+": ("+str(n)+")"
			
			#motif_a = 'NRUGAUG'
			motif_a = 'RUGAUG'.replace('U','T')
			motif_b = 'CUGA'.replace('U','T')
			
			#k=7
			#n=10
			#i[0] = 0,7
			#i[0] = 1,8
			#i[0] = 2,9
			#i[0] = 3,10
			
			k = len(motif_a)
			l = len(motif_b)
			#n = 10
			
			for i in range(n-max(k,l)+1):
				chunk1 = ref.fetch(chromosome,i,i+k)
				if match(motif_a,chunk1):
					print chromosome+":"+str(i)+"\tC-box"
				
				chunk2 = ref.fetch(chromosome,i,i+l)
				if match(motif_b,chunk2):
					print chromosome+":"+str(i)+"\tD-box"
				
			#./find.py > hg19_c_d_boxes.txt 
			
			#print
			#print match(motif_a,"A"+"TGATG"),"- should be True"
			#print match(motif_a,"G"+"TGATG"),"- should be True"
			#print match(motif_a,"C"+"TGATG"),"- should be False"
			#print match(motif_a,"T"+"TGATG"),"- should be False"
			
