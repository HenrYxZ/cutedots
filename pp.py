#!/usr/bin/env python

# Copyright 2012 Esteban Hurtado
#
# This file is part of Cutedots.
#
# Cutedots is distributed under the terms of the Reciprocal Public License 1.5.
#
# Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.
#
# You should have received a copy of the Reciprocal Public License along with
# Cutedots. If not, see <http://http://opensource.org/licenses/rpl-1.5>.



desc = """Preprocess and trajectorize 3D motion data. Performs the following tasks

    [1] Read C3D file
    [2] Find average head (solid body) pattern. The three points
        with 'H1' prefix in their label are considered to be head 1;
        'H2' for head 2.
    [3] Remove trajectorization information. New representation has
        no tracks but a list of points for each frame.
    [4] Remove duplicate points.
    [5] Track heads, labeling them track1 and track2.
    [6] Find continuous trajectories for later labeling. Final
        representation is a list of trajectories (see class 'traj' in
        file 'ngpp.pyx').
    [7] Write trajectorized data set (see class 'traj_data' in
        file 'ngpp.pyx') using input filename with '.qtd' appended.
"""

import raw_curves
import sys

import argparse
parser = argparse.ArgumentParser(
    description=desc)
parser.add_argument('-r',  action='store_true', help='Rotate 90 deg. horizontally')
parser.add_argument('file', help="Input C3D data file")
args = parser.parse_args()

raw_curves.pp(args.file, args.r)
