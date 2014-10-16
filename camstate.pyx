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

cdef float EYE_SEP = 80

# Camera position and orientation
#################################

cdef class CamState:
    def __init__(CamState self):
        self.reset()

    cdef public void reset(CamState self):
        self.px = 0
        self.py = -2000
        self.pz = 1200
        self.rx = 0
        self.ry = 0
        self.rz = 600
        self.half_eye_sep = EYE_SEP / 2.0
        self.rotation = 0.0
