#!/usr/bin/env python

"""
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2016 Youri Hoogstrate

This file is part of segmentation-fold

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

import sys,argparse,textwrap,datetime,click
from segmentation_fold_utils import __version__, __author__, __homepage__

from segmentation_fold_utils.FindBoxes import FindBoxes
from segmentation_fold_utils.ExtractBoxedSequences import ExtractBoxedSequences


@click.version_option(__version__)
@click.group()
def CLI():
	pass


@CLI.command(name='energy-estimation')
def CLI_energy_estimation():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="For more info please visit:\n<https://github.com/yhoogstrate/segmentation-fold>")
	
	parser.add_argument("-t","--temp-dir",default="/tmp",help="Directory in which the temporary files will be created")
	parser.add_argument("-s","--segmentation-fold",default="segmentation-fold",help="Path of the binary")
	parser.add_argument("-x","--xml-file",default="/usr/local/share/segmentation-fold/segments.xml",help="Location of segments.xml")
	
	parser.add_argument("-r","--randomize",default=0,type=int,help="Randomize or shuffle the sequence -r times; default 0 means this function is disabled")
	
	return parser.parse_args()


@CLI.command(name='scan-for-segments')
def CLI_scan_for_segments():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="For more info please visit:\n<https://github.com/yhoogstrate/segmentation-fold>")
	
	parser.add_argument("-t","--temp-dir",default="/tmp",help="Directory in which the temporary files will be created")
	parser.add_argument("-T","--threads",default=2,type=int,help="Used threads per spawned instance of segmentation-fold")
	parser.add_argument("-s","--segmentation-fold",default="segmentation-fold",help="Path of the binary")
	parser.add_argument("-x","--xml-file",default="/usr/local/share/segmentation-fold/segments.xml",help="Location of segments.xml")
	parser.add_argument("-p","--precision",choices=['quick','in-depth'],default="in-depth",help="Precision level 'quick' will scan quickly whether the structure changes and 'in-depth' will try to find the transition levels.")
	
	parser.add_argument("fasta_file",help="FASTA file containing the sequences that will be scanned for the presence of segments from --xml-file")
	
	return parser.parse_args()


@CLI.command(name='cd-box',short_help='Scans through a sequence for subsequences that may contain C/D-box K-turns')
@click.argument('input_fasta_file', type=click.Path(exists=True))
@click.argument('output_bed_file', type=click.File('w'))
@click.option('--box1',default='NRUGAUG',help="The first motif/box (default is C-box: NRUGAUG)")
@click.option('--box2',default='CUGA',help="The second motif/box (default is D-box: CUGA)")
@click.option('--forward/--no-foward',help="Scan in the forward strand of the genome (default: True / --forward)")
@click.option('--reverse/--no-foward',help="Scan in the reverse complement strand of the genome (default: True / --reverse)")
@click.option('--inner-dist','-d',type=int, default=250,help="The maximal distance between the boxes (default=250).")
def CLI_scan_for_cd_box_kturns(input_fasta_file,box1,box2,forward,reverse,output_bed_file):
	boxes = FindBoxes(input_fasta_file,box1,box2,forward,reverse,inner_dist)
	boxes.run(output_bed_file)
	#sequences = ExtractBoxedSequences(fasta_input_file,bed_input_file,fasta_output_file,max_inner_dist,bp_extension)
	#sequences.run(fasta_output_file)


@CLI.command(name='find-boxes',short_help='find all occurances of two given boxes (sequence motifs) within a FASTA file')
@click.argument('fasta_input_file',  type=click.Path(exists=True))
@click.argument('bed_output_file', type=click.File('w'))
@click.option('--box1','-c',default='NRUGAUG',help="Sequence of box1 (default = C-box: 'NRUGAUG')")
@click.option('--box2','-d',default='CUGA',help="Sequence of box2 (default = D-box: 'CUGA')")
@click.option('--forward/--no-forward',default=True,help="Search in the forward direction of the reference sequence")
@click.option('--reverse/--no-reverse',default=True,help="Search in the reverse complement of the reference sequence")
def CLI_extract_boxed_sequences(fasta_input_file,bed_output_file,box1,box2,forward,reverse):
    boxes = FindBoxes(fasta_input_file,box1,box2,forward,reverse)
    boxes.run(bed_output_file)


@CLI.command(name='extract-boxed-sequences',short_help='bed_input_file has to be  created with \'find-box\' as part of this utility')
@click.argument('fasta_input_file',  type=click.Path(exists=True))
@click.argument('bed_input_file', type=click.File('r'))
@click.argument('fasta_output_file', type=click.File('w'))
@click.option('--max-inner-dist','-d',type=int, default=250,help="Maximal distance between the boxes (default=250bp)")
@click.option('--bp-extension','-e',type=int, default=10,help="Extend the extraced boxed sequences with this umber of bases (default: 10bp)")
def CLI_extract_boxed_sequences(fasta_input_file,bed_input_file,fasta_output_file,max_inner_dist,bp_extension):
    sequences = ExtractBoxedSequences(fasta_input_file,bed_input_file,fasta_output_file,max_inner_dist,bp_extension)
    sequences.run(fasta_output_file)

