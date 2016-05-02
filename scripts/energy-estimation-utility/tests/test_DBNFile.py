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

from segmentation_fold_utils.DBNFile import DBNFile

def get_n_lines(filename):
    i = 0
    with open(filename,"r") as fh:
        for line in fh:
            i += 1
    return i

class TestDBNFile(unittest.TestCase):
    def test_01(self):
        dbn_input_file = 'DBNFile.test_01.in.dbn'
        bam_input_file = 'DBNFile.test_01.in.bam'
        dbn_output_file = 'DBNFile.test_01.out.dbn'
        
        i = 0
        structures = DBNFile(open("tests/test-data/"+dbn_input_file,"r"),True)# open as energy-estimation dbn
        
        for sequence in structures.parse():
            if i == 0:
                self.assertEqual(sequence['name'],">chr1:10-21 x unknown-01")
                self.assertEqual(sequence['sequence'],"GGGGAAACCCC")
                self.assertEqual(len(sequence['transitions']),2)
            elif i == 1:
                self.assertEqual(sequence['name'],">chr1:25-36 x unknown-01")
                self.assertEqual(sequence['sequence'],"AAAAAAAAAAA")
                self.assertEqual(len(sequence['transitions']),0)
            elif i == 2:
                self.assertEqual(sequence['name'],">chr1:45-56 x unknown-01")
                self.assertEqual(sequence['sequence'],"AAAAAAAAAAA")
                self.assertEqual(len(sequence['transitions']),0)
            
            i += 1
        
        self.assertEqual(i,3)
        
        structures.annotate_read_counts('tests/test-data/'+bam_input_file,'>.*?(chr[^:]):([0-9]+)-([0-9]+)',open(dbn_output_file,"w"))
        
        self.assertTrue(filecmp.cmp(dbn_output_file,"tests/test-data/"+dbn_output_file))
    
    def test_02(self):
        dbn_input_file = 'DBNFile.test_02.in.dbn'
        bed_input_file = 'DBNFile.test_02.in.bed'
        dbn_output_file_o = 'DBNFile.test_02.out.o.dbn'
        dbn_output_file_n = 'DBNFile.test_02.out.n.dbn'
        
        structures = DBNFile(open("tests/test-data/"+dbn_input_file,"r"),True)# open as energy-estimation dbn
        structures.filter_annotated_entries( \
            '>.*?(chr[^:]):([0-9]+)-([0-9]+)', \
            open("tests/test-data/"+bed_input_file,"r"), \
            open(dbn_output_file_o,"w"), \
            open(dbn_output_file_n,"w"), \
        )
        
        self.assertTrue(filecmp.cmp(dbn_output_file_o,"tests/test-data/"+dbn_output_file_o))
        self.assertTrue(filecmp.cmp(dbn_output_file_n,"tests/test-data/"+dbn_output_file_n))

def main():
    unittest.main()

if __name__ == '__main__':
    main()
