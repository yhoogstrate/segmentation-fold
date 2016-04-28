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

class TestXMLFile(unittest.TestCase):
	def test_01(self):
		input_xml_file = 'XMLFile.test_01.in.xml'
		input_fasta_file = None
		
		xml = XMLFile("tests/test-data/"+input_xml_file)
		for a,b,c,d in xml.get_combinations(input_fasta_file):
			print "-----------------------------------"
			print str(a)
			print str(b)
			print str(c)
			print str(d)
			print "-----------------------------------\n\n\n\n"
	
	def test_02(self):
		input_xml_file = 'XMLFile.test_01.in.xml'
		input_fasta_file = 'XMLFile.test_02.in.fa'
		
		xml = XMLFile("tests/test-data/"+input_xml_file)
		for a,b,c,d in xml.get_combinations("tests/test-data/"+input_fasta_file):
			print "-----------------------------------"
			print str(a)
			print str(b)
			print str(c)
			print str(d)
			print "-----------------------------------\n\n\n\n"

def main():
	unittest.main()

if __name__ == '__main__':
	main()
