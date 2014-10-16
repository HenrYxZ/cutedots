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
# Cutedots. If not, see <http://opensource.org/licenses/rpl-1.5>.

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

if os.name == 'nt':
    display_libs=['opengl32','glu32']
else:
    display_libs=['GL', 'GLU']
    
setup(cmdclass = {'build_ext': build_ext}, ext_modules = [
        Extension("dotsdisplay",
                  ['dotsdisplay.pyx', 'dotsdisplay.pxd', 'glextern.pxd', 'dots_drawing.c'],
                  libraries=display_libs),
        Extension("camstate", ["camstate.pyx", "camstate.pxd"]),
        Extension("modelstate", ["modelstate.pyx", "modelstate.pxd"]),
        Extension("c3dformat", ["c3dformat.pyx"]),
        Extension("trajdata", ["trajdata.pyx"]),
        Extension("dotsio", ["dotsio.pyx"]),
        Extension("modelops", ["modelops.pyx"]),
        Extension("preprocess", ["preprocess.pyx"]),
        ])
