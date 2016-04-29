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

from segmentation_fold_utils.XMLFile import XMLFile

def get_n_lines(filename):
    i = 0
    with open(filename,"r") as fh:
        for line in fh:
            i += 1
    return i

class TestXMLFile(unittest.TestCase):
    def test_01(self):
        input_xml_file = 'XMLFile.test_01.in.xml'
        input_fasta_file = None
        dbn_output_file = 'XMLFile.test_01.out.dbn'
        
        i = 0
        xml = XMLFile(open("tests/test-data/"+input_xml_file,"r"))
        for sequence_name, sequence, segment_name, segment in xml.get_combinations(input_fasta_file):
            if i == 0:
                self.assertEqual(sequence_name,"C/D-box snoRNA")
                self.assertEqual(sequence,"GCUCUGACCGAAAGGCGUGAUGAGC")
                self.assertEqual(segment_name,"Kt-CD-box.CGU")
                self.assertEqual(segment['5'],"CGUGAU")
                self.assertEqual(segment['x'],'   :::')
                self.assertEqual(segment['3'],'   AGU')
            elif i == 1:
                self.assertEqual(sequence_name,"Artificial double C/D K-turn construct")
                self.assertEqual(sequence,"GGGAGUCUUGUGAUGAGAAGUACUGGAUCUGAAGUAGCCCUUUUUGGGCUACUUGUGAUGAAACACUCAUGGUCUGAAGACUCCC")
                self.assertEqual(segment_name,"Kt-CD-box.UGU")
                self.assertEqual(segment['5'],"UGUGAU")
                self.assertEqual(segment['x'],'   :::')
                self.assertEqual(segment['3'],'   AGU')
            i += 1
        
        self.assertEqual(i,2)
        
        xml.estimate_energy("/tmp","segmentation-fold",1,0,0,input_fasta_file,open(dbn_output_file,"w"))
        self.assertEqual(get_n_lines(dbn_output_file),7)
    
    def test_02(self):
        input_xml_file = 'XMLFile.test_01.in.xml'
        input_fasta_file = "tests/test-data/"+'XMLFile.test_02.in.fa'
        dbn_output_file = 'XMLFile.test_02.out.dbn'
        
        i = 0
        xml = XMLFile(open("tests/test-data/"+input_xml_file,"r"))
        for sequence_name, sequence, segment_name, segment in xml.get_combinations(input_fasta_file):
            if i == 0:
                self.assertEqual(sequence_name,"C/D-box snoRNA")
                self.assertEqual(sequence,"GCUCUGACCGAAAGGCGUGAUGAGC")
                self.assertEqual(segment_name,"Kt-CD-box.CGU")
                self.assertEqual(segment['5'],"CGUGAU")
                self.assertEqual(segment['x'],'   :::')
                self.assertEqual(segment['3'],'   AGU')
            if i == 1:
                self.assertEqual(sequence_name,"C/D-box snoRNA")
                self.assertEqual(sequence,"GCUCUGACCGAAAGGCGUGAUGAGC")
                self.assertEqual(segment_name,"Kt-CD-box.UGU")
                self.assertEqual(segment['5'],"UGUGAU")
                self.assertEqual(segment['x'],'   :::')
                self.assertEqual(segment['3'],'   AGU')
            elif i == 2:
                self.assertEqual(sequence_name,"Artificial double C/D K-turn construct")
                self.assertEqual(sequence,"GGGAGUCUUGUGAUGAGAAGUACUGGAUCUGAAGUAGCCCUUUUUGGGCUACUUGUGAUGAAACACUCAUGGUCUGAAGACUCCC")
                self.assertEqual(segment_name,"Kt-CD-box.CGU")
                self.assertEqual(segment['5'],"CGUGAU")
                self.assertEqual(segment['x'],'   :::')
                self.assertEqual(segment['3'],'   AGU')
            elif i == 3:
                self.assertEqual(sequence_name,"Artificial double C/D K-turn construct")
                self.assertEqual(sequence,"GGGAGUCUUGUGAUGAGAAGUACUGGAUCUGAAGUAGCCCUUUUUGGGCUACUUGUGAUGAAACACUCAUGGUCUGAAGACUCCC")
                self.assertEqual(segment_name,"Kt-CD-box.UGU")
                self.assertEqual(segment['5'],"UGUGAU")
                self.assertEqual(segment['x'],'   :::')
                self.assertEqual(segment['3'],'   AGU')
            i += 1
        
        self.assertEqual(i,4)
        
        xml.estimate_energy("/tmp","segmentation-fold",1,0,0,input_fasta_file,open(dbn_output_file,"w"))
        self.assertEqual(get_n_lines(dbn_output_file),11)
    
    def test_03(self):
        input_xml_file = 'XMLFile.test_01.in.xml'
        input_fasta_file = None
        dbn_output_file = 'XMLFile.test_03.out.dbn'
        
        i = 0
        xml = XMLFile(open("tests/test-data/"+input_xml_file,"r"))
        xml.estimate_energy("/tmp","segmentation-fold",1,0,50,input_fasta_file,open(dbn_output_file,"w"))
        self.assertTrue(get_n_lines(dbn_output_file) >= 200)
    

def main():
    unittest.main()

if __name__ == '__main__':
    main()
