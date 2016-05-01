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
from segmentation_fold_utils.XMLFile import XMLFile


@click.version_option(__version__)
@click.group()
def CLI():
    pass


@CLI.command(name="estimate-energy",short_help="Estimates whether a certain Segment(Loop) is present and for which delta-G this transistion takes place")
@click.option("--temp-dir","-t",type=click.Path(exists=True),default="/tmp",help="Ppath in which temporary conifig files will be stored (default: /tmp)")
@click.option("--segmentation-fold","-s",default="segmentation-fold",help="Location of segmentatio-fold binary (default: segmentation-fold)")
@click.option("--xml-file","-x", type=click.File('r'),default="/usr/local/share/segmentation-fold/segments.xml",help="Location of segments.xml (default: /usr/local/share/segmentation-fold/segments.xml)")
@click.option("--threads","-T",default=2,type=int,help="Number of threads per spawned instance of segmentation-fold")
@click.option("--precision","-p",default=0.05,help="Minimal difference for binary split - the smaller this value the slower. if this value equals 0, the difference is set to infinity (default: 0.05)")
@click.option("--randomize","-r",default=0,type=int,help="Shuffle each sequence this many times and predict energy of shuffled sequence(s) (default: 0, 0 means disabled)")
@click.option("--sequences-from-fasta-file","-f", type=click.Path(exists=True),default=None,help="Use sequences from a FASTA file instead of using the XML file that also contains the segments. In XML files you can explicitly link one Segment(Loop) to one particular sequence instead of doing n*n comparisons (default: None)")
@click.argument('dbn_output_file', type=click.File('w'))
def CLI_energy_estimation(temp_dir,segmentation_fold,xml_file,threads,precision,randomize,sequences_from_fasta_file,dbn_output_file):
    xml = XMLFile(xml_file)
    xml.estimate_energy(temp_dir,segmentation_fold,threads,precision,randomize,sequences_from_fasta_file,dbn_output_file)


@CLI.command(name='find-boxes',short_help='Finds all occurances of two given boxes (sequence motifs) within a FASTA file')
@click.argument('fasta_input_file',  type=click.Path(exists=True))
@click.argument('bed_output_file', type=click.File('w'))
@click.option('--box1','-c',default='NRUGAUG',help="Sequence of box1 (default = C-box: 'NRUGAUG')")
@click.option('--box2','-d',default='CUGA',help="Sequence of box2 (default = D-box: 'CUGA')")
@click.option('--forward/--no-forward',default=True,help="Search in the forward direction of the reference sequence")
@click.option('--reverse/--no-reverse',default=True,help="Search in the reverse complement of the reference sequence")
def CLI_extract_boxed_sequences(fasta_input_file,bed_output_file,box1,box2,forward,reverse):
    boxes = FindBoxes(fasta_input_file,box1,box2,forward,reverse)
    boxes.run(bed_output_file)


@CLI.command(name='extract-boxed-sequences',short_help='Extracts boxes sequences from bed_input_file which has to be created with \'find-box\', part of this utility')
@click.argument('fasta_input_file',  type=click.Path(exists=True))
@click.argument('bed_input_file', type=click.File('r'))
@click.argument('fasta_output_file', type=click.File('w'))
@click.option('--max-inner-dist','-d',type=int, default=250,help="Maximal distance between the boxes (default=250bp)")
@click.option('--bp-extension','-e',type=int, default=10,help="Extend the extraced boxed sequences with this umber of bases (default: 10bp)")
def CLI_extract_boxed_sequences(fasta_input_file,bed_input_file,fasta_output_file,max_inner_dist,bp_extension):
    sequences = ExtractBoxedSequences(fasta_input_file,bed_input_file,fasta_output_file,max_inner_dist,bp_extension)
    sequences.run(fasta_output_file)


@CLI.command(name='add-read-counts',short_help='Add read counts falling within a region contained in the fasta header of the dbn file')
@click.argument('dbn_input_file', type=click.Path(exists=True))
@click.argument('bam_input_file', type=click.Path(exists=True))
@click.argument('dbn_output_file', type=click.File('w'))
@click.option('--regex','-r', default=">.*?(chr[^:]):([0-9]+)-([0-9]+)",help="Regex to capture the targeted location")
# possible more sam/bam flag requirements, e.g. multimap limits etc.
def CLI_add_read_counts(fasta_input_file,bed_input_file,fasta_output_file,max_inner_dist,bp_extension):
    sequences = ExtractBoxedSequences(fasta_input_file,bed_input_file,fasta_output_file,max_inner_dist,bp_extension)
    sequences.run(fasta_output_file)

