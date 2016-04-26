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


import unittest,logging,sys,filecmp
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from segmentation_fold_utils.ExtractBoxedSequences import ExtractBoxedSequences

class TestExtractBoxedSequences(unittest.TestCase):
    def test_01(self):
        fasta_input_file = 'ExtractBoxedSequences.test_01.in.fa'
        bed_input_file =   "ExtractBoxedSequences.test_01.in.bed"
        fasta_output_file = "ExtractBoxedSequences.test_01.out.fa"

        sequences = ExtractBoxedSequences('tests/test-data/'+fasta_input_file,open('tests/test-data/'+bed_input_file,"r"),open('tests/test-data/'+fasta_output_file,'r'),90,0)
        sequences.run(open(fasta_output_file,"w"))
        
        self.assertTrue(filecmp.cmp(fasta_output_file,"tests/test-data/"+fasta_output_file) )

    def test_02(self):
        fasta_input_file = 'ExtractBoxedSequences.test_02.in.fa'
        bed_input_file =   "ExtractBoxedSequences.test_02.in.bed"
        fasta_output_file = "ExtractBoxedSequences.test_02.out.fa"

        sequences = ExtractBoxedSequences('tests/test-data/'+fasta_input_file,open('tests/test-data/'+bed_input_file,"r"),open('tests/test-data/'+fasta_output_file,'r'),90,0)
        sequences.run(open(fasta_output_file,"w"))
        
        self.assertTrue(filecmp.cmp(fasta_output_file,"tests/test-data/"+fasta_output_file) )



def main():
    unittest.main()

if __name__ == '__main__':
    main()

