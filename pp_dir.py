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

import os
import sys
from subprocess import call

base_dir = sys.argv[1]
files = os.listdir(base_dir)

for fn in files:
    if fn[-4:] == ".c3d":
        call([ "python", "pp.py", os.path.join(base_dir, fn) ])
#        call([ "python", "join.py", os.path.join(base_dir, fn+".pickle") ])
