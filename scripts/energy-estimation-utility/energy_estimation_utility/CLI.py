#!/usr/bin/env python

"""
@file scripts/energy-estimation-utility/energy_estimation_utility/CLI.py

@date 2015-07-29

@author Youri Hoogstrate

@section LICENSE
<PRE>
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2015 Youri Hoogstrate

This file is part of segmentation-fold and originally taken from
yh-kt-fold.

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

import sys,argparse,textwrap,datetime

def CLI():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="For more info please visit:\n<https://github.com/yhoogstrate/segmentation-fold>")
	
	parser.add_argument("-t","--temp-dir",default="/tmp",help="Directory in which the temporary files will be created")
	parser.add_argument("-s","--segmentation-fold",default="segmentation-fold",help="Path of the binary")
	parser.add_argument("-x","--xml-file",default="/usr/local/share/segmentation-fold/segments.xml",help="Location of segments.xml")
	
	return parser.parse_args()



def CLI_scan_for_segments():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="For more info please visit:\n<https://github.com/yhoogstrate/segmentation-fold>")
	
	parser.add_argument("-t","--temp-dir",default="/tmp",help="Directory in which the temporary files will be created")
	parser.add_argument("-s","--segmentation-fold",default="segmentation-fold",help="Path of the binary")
	parser.add_argument("-x","--xml-file",default="/usr/local/share/segmentation-fold/segments.xml",help="Location of segments.xml")
	parser.add_argument("-p","--precision",choices=['quick','in-depth'],default="quick",help="Precision level 'quick' will scan quickly whether the structure changes and 'in-depth' will try to find the transition levels.")
	
	parser.add_argument("fasta_file",help="FASTA file containing the sequences that will be scanned for the presence of segments from --xml-file")
	
	return parser.parse_args()