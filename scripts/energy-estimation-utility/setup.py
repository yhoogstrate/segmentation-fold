#!/usr/bin/env python

"""
@file scripts/energy-estimation/setup.py

@author Youri Hoogstrate

@section LICENSE
<PRE>
segmentation-fold can predict RNA 2D structures including K-turns.
Copyright (C) 2012-2016 Youri Hoogstrate

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

import segmentation_fold_utils

#from setuptools.command.install import install
from distutils.core import setup
from setuptools import setup, find_packages

import os

#class RunMarkdown(install):
#	def run(self):
#		install.run(self)
#		os.system('markdown bin/README.md > bin/README.html 2> /dev/null')

setup(name='segmentation-fold-utils',
		version=segmentation_fold_utils.__version__,
		description='energy estimation utils for segmentation-fold',
		author=segmentation_fold_utils.__author__,
		maintainer=segmentation_fold_utils.__author__,
		license=segmentation_fold_utils.__license__,
		url=segmentation_fold_utils.__homepage__,
		scripts=["bin/segmentation-fold-utils"],
		packages=['segmentation_fold_utils'],
		test_suite="tests",
		install_requires=[
			'HTSeq >= 0.6.1',
			'pysam >= 0.8.0'
		],
		classifiers=[
			'Environment :: Console',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
			'Operating System :: OS Independent'
			'Topic :: Scientific/Engineering',
			'Topic :: Scientific/Engineering :: Bio-Informatics',
			],
#		cmdclass={'install':RunMarkdown}
	)
