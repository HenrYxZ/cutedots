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

cdef class CamState:
    cdef public float px        # Camera position
    cdef public float py
    cdef public float pz
    cdef public float rx       # Camara reference point
    cdef public float ry
    cdef public float rz
    cdef public float half_eye_sep    # Half of eye separation
    cdef public float rotation

    cdef public void reset(CamState self)

