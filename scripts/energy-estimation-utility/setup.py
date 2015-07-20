#!/usr/bin/env python

"""
@file scripts/energy-estimation/bin/energy-estimation

@date 2015-07-15

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

import energy_estimation_utility

from distutils.core import setup
from setuptools import setup, find_packages

setup(name='energy-estimation-utility',
		version=energy_estimation_utility.__version__,
		description='energy-estimation-utility for segmentation-fold',
		author=energy_estimation_utility.__author__,
		maintainer=energy_estimation_utility.__author__,
		license=energy_estimation_utility.__license__,
		url=energy_estimation_utility.__homepage__,
		scripts=["bin/energy-estimation-utility"],
		packages=['energy_estimation_utility'],
		classifiers=[
			'Environment :: Console',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
			'Operating System :: OS Independent'
			'Topic :: Scientific/Engineering',
			'Topic :: Scientific/Engineering :: Bio-Informatics',
			],
	)
