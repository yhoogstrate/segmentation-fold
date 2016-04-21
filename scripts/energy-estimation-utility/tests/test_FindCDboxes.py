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


import unittest,logging,sys
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from segmentation_fold_utils.FindCDboxes import FindCDboxes

class TestFindCDboxes(unittest.TestCase):
	def test_01(self):
		output_file = "TestFindCDboxes.test_01.txt"
		boxes = FindCDboxes('tests/test-data/FindCDboxes.genome.txt','NRUGAUG','CUGA',True,True,250)
		boxes.run(output_file)

def main():
	unittest.main()

if __name__ == '__main__':
	main()
