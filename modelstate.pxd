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

from camstate cimport *
from trajdata import TrajData

# Display modes
###############

cdef enum display_mode:
     display_frame, display_all, display_curve

cdef class ModelState:
    cdef public object data           # Reference to data
    cdef public int frame             # current frame
    cdef public int traj              # current trajectory
    cdef public CamState cam          # Camera
    cdef public int anaglyph3d        # Whether to render anaglyph 3D
    cdef public float bodyTrans       # Body transparency
    cdef public int width             # Screen dimensions
    cdef public int height
    cdef public float aspect
    cdef public list trash            # Stores deleted curves

    # For skeleton lines
    cdef public dict trajMap


