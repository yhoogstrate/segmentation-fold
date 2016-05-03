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


import os,logging
import pysam


class FindBoxes:
	"""detecting C boxes in chr22 takes ~7min
	
	results are in BED format: "0-based, half open" - i.e. first position is: chr1\t0\t1
	"""
	
	logger = logging.getLogger("segmentation-fold-utils::FindBoxes")
	
	def __init__(self,genome,box1,box2,search_fwd,search_rev):
		self.genome = genome
		self.box1 = box1.upper().replace('U','T')#.strip('N')
		self.box2 = box2.upper().replace('U','T')#.strip('N')
		self.k = len(self.box1)
		self.l = len(self.box2)
		self.m = min(self.k, self.l)
		self.search_fwd = search_fwd
		self.search_rev = search_rev
		
		if self.search_rev:
			self.box1r = self.reverse_complement(self.box1)
			self.box2r = self.reverse_complement(self.box2)
	
	def reverse_complement(self, sequence):
		"""Instead of RC' the entire genome, as what would biologically
		make sense because of the double stranded DNA, we 'reverse complement'
		the box motif.
		
		------ fwd ------>
		aaaaaaaCTGAaaaaaaa
		       ||||
		motif: CTGA
		
		<------ rev ------
		tttttttTCAGttttttt
		       ||||
		motif: TCAG
		
		---------------------------------
		motif containing R (A|G)
		
		------ fwd ------>
		aaaaaaaCAGAaaaaaaa <- match1
		aaaaaaaCGGAaaaaaaa <- match1
		       ||||
		motif: CRGA
		
		<------ rev ------
		tttttttTC{rc:R}Gttttttt
		       ||||
		motif: TCAG
		
		(reverse searching should also start with the highest position)
		
		
		"""
		rc = ''
		tt = {
			'N':'N',
			'A':'T',
			'C':'G',
			'T':'A',
			'G':'C',
			
			'R':'Y',# R = (A||G), rc(A||G) = (T||C) = Y
			'Y':'R',
			'K':'M',
			'M':'K',# M = (A||C), rc(A||C) = (T||G) = K
			'S':'S',# S = (G||C), rc(G||C) = (C||G) = S
			'W':'W',# W = (A||U), rc(A||U) = (U||A) = W
			
			'B':'V',# B = (C||G||T), rc(C||G||T) = (G||C||A) = V
			'D':'H',# D = AGT, rc(AGT) = TCA = H
			'H':'D',# H = ACT, rc(ACT) = TGA = D
			'V':'B'# V = ACG, rc(ACG) = TGC = B
		}
		
		for base in sequence[::-1]:
			rc += tt[base]
		
		return rc
	
	def match_base(self,base_pattern, base_query):
		if base_query == "A":
			return base_pattern in "ARMWDHVN"
		elif base_query == "G":
			return base_pattern in "GRKSBDVN"
		elif base_query == "C":
			return base_pattern in "CYMSBHVN"
		elif base_query == "T":
			return base_pattern in "TYKWBDHN"
		else:
			return False
	
	def match(self,pattern, query):
		if len(pattern) != len(query):
			return False
		
		for i in range(len(pattern)):
			if not self.match_base(pattern[i], query[i]):
				return False
		return True
	
	def check_faid_out_of_date(self,output_file):
		faid_index = output_file+".fai"
		if os.path.isfile(faid_index) and os.path.isfile(output_file):
			if os.path.getmtime(faid_index) < os.path.getmtime(output_file):
				self.logger.info("The faid index of "+output_file+" is older than the file itself. Removing the index.")
				os.remove(faid_index)
	
	def run(self,fh):
		self.check_faid_out_of_date(self.genome)
		ref = pysam.FastaFile(self.genome)
		for chromosome in ref.references:
			# Look fwd
			n = ref.get_reference_length(chromosome)
			self.logger.debug("scanning "+chromosome+" ("+str(n)+" bases)")
			
			for i in range(n-self.m+1):
				k = i+self.k
				l = i+self.l
				
				chunk1 = ref.fetch(chromosome,i,k).upper().replace('U','T')
				chunk2 = ref.fetch(chromosome,i,l).upper().replace('U','T')
				
				if self.search_fwd:
					if self.match(self.box1,chunk1):
						fh.write(chromosome+"\t"+str(i)+"\t"+str(k)+"\tbox1-f:"+self.box1+"\t0\t+\n")
					if self.match(self.box2,chunk2):
						fh.write(chromosome+"\t"+str(i)+"\t"+str(l)+"\tbox2-f:"+self.box2+"\t0\t+\n")
				
				if self.search_rev:
					if self.match(self.box1r,chunk1):
						fh.write(chromosome+"\t"+str(i)+"\t"+str(k)+"\tbox1-r:"+self.box1r+"\t0\t-\n")
					if self.match(self.box2r,chunk2):
						fh.write(chromosome+"\t"+str(i)+"\t"+str(l)+"\tbox2-r:"+self.box2r+"\t0\t-\n")
		
		fh.close()
