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


class Test_filter_annotated_entries(unittest.TestCase):
    def test_01(self):
        """
        Do the unit test via the command line
         - Requires segmentation-fold-utils to be installed
        """
        dbn_input_file = 'DBNFile.test_02.in.dbn'
        bed_input_file = 'DBNFile.test_02.in.bed'
        dbn_output_file_o = 'DBNFile.test_02.out.o.dbn'
        dbn_output_file_n = 'DBNFile.test_02.out.n.dbn'
        
        command = ["segmentation-fold-utils", \
                    "filter-annotated-entries", \
                    "--regex",">.*?(chr[^:]):([0-9]+)-([0-9]+)",
                    "tests/test-data/"+dbn_input_file, \
                    "tests/test-data/"+bed_input_file, \
                    dbn_output_file_o, \
                    dbn_output_file_n]
        
        self.assertEqual(subprocess.call(command) , 0)
        
        self.assertTrue(filecmp.cmp(dbn_output_file_o,"tests/test-data/"+dbn_output_file_o))
        self.assertTrue(filecmp.cmp(dbn_output_file_n,"tests/test-data/"+dbn_output_file_n))


def main():
    unittest.main()

if __name__ == '__main__':
    main()
