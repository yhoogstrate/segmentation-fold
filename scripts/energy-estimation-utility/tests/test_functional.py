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


import unittest,logging,sys,filecmp,subprocess
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)


class Test_find_boxes(unittest.TestCase):
    def test_01(self):
        """
        Do the unit test via the command line
         - Requires segmentation-fold-utils to be installed
        """
        input_file = "functional.test_01.in.fa"
        output_file = "functional.test_01.tmp.bed"
        
        command = ["segmentation-fold-utils", \
                    "find-boxes", \
                    "--box1","NRUGAUG",
                    "--box2","CUGA",
                    "--forward",
                    "--reverse",
                    "tests/test-data/"+input_file, \
                    output_file]
        
        self.assertEqual(subprocess.call(command) , 0)
        #self.assertTrue(filecmp.cmp(output_file,"tests/test-data/"+output_file))
        
        fasta_input_file = input_file
        bed_input_file = output_file
        fasta_output_file = "functional.test_01.tmp.fa"
        
        command = ["segmentation-fold-utils",
                    "extract-boxed-sequences",
                    "--max-inner-dist","90",
                    "--bp-extension","3",
                    "tests/test-data/"+fasta_input_file,
                    bed_input_file,
                    fasta_output_file]
        
        self.assertEqual(subprocess.call(command) , 0)
        #self.assertTrue(filecmp.cmp(fasta_output_file,"tests/test-data/"+fasta_output_file))
        
        input_xml_file = 'XMLFile.test_01.in.xml'
        input_fasta_file = fasta_output_file
        dbn_output_file = 'functional.test_01.tmp.dbn'
        command = ["segmentation-fold-utils",
                   "estimate-energy",
                   "--xml-file","tests/test-data/"+input_xml_file,
                   "--threads","1",
                   "--precision","0",
                   "--randomize","0",
                   "--sequences-from-fasta-file",input_fasta_file,
                   dbn_output_file]
        
        self.assertEqual(subprocess.call(command) , 0)
        ##self.assertEqual(get_n_lines(dbn_output_file),7)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
