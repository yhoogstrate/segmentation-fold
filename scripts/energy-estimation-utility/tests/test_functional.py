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


    def test_02(self):
        """
        Do the unit test via the command line
         - Requires segmentation-fold-utils to be installed
        """
        input_file_fa = "workflow-test_cd-box_kturns.fa"
        input_file_xml = "workflow-test_cd-box_kturns.xml"
        output_file_f1 = "workflow-test_f1.dbn"
        output_file_f2 = "workflow-test_f2.dbn"
        
        command_01 = ['segmentation-fold-utils',
                      'fix-fasta-headers',
                      'tests/test-data/'+input_file_fa,
                      'fasta_fixed_headers.fa']
        command_02 = ['segmentation-fold-utils',
                      'find-boxes',
                      '--box1',
                      'NRUGAUG',
                      '--box2',
                      'CUGA',
                      '--forward',
                      '--reverse',
                      'fasta_fixed_headers.fa',
                      'detected_boxes.bed']
        command_03 = ['segmentation-fold-utils',
                      'extract-boxed-sequences',
                      '--max-inner-dist',
                      '110',
                      '--bp-extension',
                      '10',
                      'fasta_fixed_headers.fa',
                      'detected_boxes.bed',
                      'sequences_containing_boxes_within_110bp.fa']
        command_04 = ['segmentation-fold-utils',
                      'estimate-energy',
                      '-T',
                      '1',
                      '-x',
                      'tests/test-data/'+input_file_xml,
                      '-p',
                      '0.1',
                      '-r',
                      '0',
                      '--sequences-from-fasta-file',
                      'sequences_containing_boxes_within_110bp.fa',
                      'energy_estimation_raw_results.dbn']
        command_05 = ['segmentation-fold-utils',
                      'filter-by-energy',
                      '--energy',
                      '-15.0',
                      'energy_estimation_raw_results.dbn',
                      output_file_f1,
                      output_file_f2]
        
        self.assertEqual(subprocess.call(command_01) , 0)
        self.assertEqual(subprocess.call(command_02) , 0)
        self.assertEqual(subprocess.call(command_03) , 0)
        self.assertEqual(subprocess.call(command_04) , 0)
        self.assertEqual(subprocess.call(command_05) , 0)
        
        
        self.assertTrue(filecmp.cmp(output_file_f1, "tests/test-data/"+output_file_f1))
        self.assertTrue(filecmp.cmp(output_file_f2, "tests/test-data/"+output_file_f2))

def main():
    unittest.main()

if __name__ == '__main__':
    main()
