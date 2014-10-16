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

import c3dformat
from numpy import size
import sys
from numpy import array
from scipy import polyfit, polyval
from sys import float_info
import os

# Use C stdlib functions for memory allocation
cdef extern from "stdlib.h":
    void free(void* ptr)
    void* malloc(size_t size)
    void* realloc (void* ptr, size_t size)

def log(msg, lf=True):
    "Log a message"
    print(msg)
    sys.stdout.flush()

    
### RAW data storage
####################

cdef class RawFrame:
    "Non-trajectorized frame"

    cdef float *data
    cdef public int numPoints
    cdef public set unusedIndices

    def __cinit__(self, arr):
        # Find number of valid points.
        self.numPoints = sum(arr[:,3] >= 0.0)
        cdef int N = self.numPoints
       # Allocate space for them.
        self.data = <float*>malloc(3*N*sizeof(float))
        # Copy them
        cdef int idx = 0
        cdef int i
        cdef int Ntotal = arr.shape[0]
        for i from 0 <= i < Ntotal:
            x,y,z,mask = arr[i,:]
            if mask >= 0.0:
                self.data[idx] = x
                idx += 1
                self.data[idx] = y
                idx += 1
                self.data[idx] = z
                idx += 1
        self.unusedIndices = set(range(self.numPoints))

    def getSinglePoint(RawFrame self, index):
        cdef float* point = self.getPoint(index)
        return (point[0], point[1], point[2])

    cdef float* getPoint(self, i):
        return &self.data[3*i]

    def joinClosePoints(RawFrame self, float max_dist=100):
        # Allocate space to mark unique points
        cdef int *uniquePoint = <int*>malloc(self.numPoints * sizeof(int))
        cdef int i,j
        for i from 0 <= i < self.numPoints: # Initially all unique
            uniquePoint[i] = 1
        # Mark according to pairwise distances
        cdef float *a, *b
        cdef float distance
        cdef int newNumPoints = self.numPoints
        for i from 0 <= i < self.numPoints:
            for j from (i+1) <= j < self.numPoints:
                a = self.getPoint(i)
                b = self.getPoint(j)
                sq_distance = \
                    (b[0]-a[0])*(b[0]-a[0]) + \
                    (b[1]-a[1])*(b[1]-a[1]) + \
                    (b[2]-a[2])*(b[2]-a[2])
                if sq_distance < max_dist:
                    if uniquePoint[j]:
                        newNumPoints -= 1                      
                    uniquePoint[j] = 0
        # Rebuild frame data
        cdef float *newData = <float*>malloc(3*newNumPoints*sizeof(float))
        cdef float *index = newData
        for i from 0 <= i < self.numPoints:
            if uniquePoint[i]:
                index[0] = self.data[3*i]
                index[1] = self.data[3*i+1]
                index[2] = self.data[3*i+2]
                index = &index[3]
        # Assign new data
        free(self.data)
        self.data = newData
        self.numPoints = newNumPoints
        self.unusedIndices = set(range(self.numPoints))
        # Free unique point marks
        free(uniquePoint)

    def __dealloc__(self):
        free(self.data)
        self.data = NULL

class RawDataReadError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return 'Raw data read error: ' + self.msg

class RawData:
    "Array of raw frames"

    @property
    def numFrames(self):
        return len(self.frames)

    def __init__(self):
        self.frames = []
        self.filename = None
        self.frameRate = 100.0

    def joinClosePoints(self, progress):
        cdef int i = 0
        for frame in self.frames:
            if i % 1000 == 0:
                progress.setValue( int(100.0*i / self.numFrames) )
                if progress.wasCanceled():
                    return                
            frame.joinClosePoints()
            i += 1
        progress.setValue(100)


### Trajectories
################

class Traj:
    "Same point at different frames."

    @property
    def endFrame(self):
        return self.beginFrame + len(self.pointData)

    @property
    def numFrames(self):
        return len(self.pointData)

    def __init__(self, int beginFrame, str name):
        self.name = name
        self.pointData = []
        self.beginFrame = beginFrame

    def addPoint(self, point):
        self.pointData.append(point)
        self.endFrame += 1

    def average(self):
        cdef float sx=0.0, sy=0.0, sz=0.0
        cdef int n = self.numFrames
        for x, y, z in self.pointData:
            sx += x;  sy += y; sz += z
        return sx/n, sy/n, sz/n

    def averageX(self):
        cdef float sx=0.0
        for p in self.pointData:
            sx += p[0]
        return sx / self.numFrames

    def hasFrame(self, int framenum):
        """Returns True if given frame is included in trajectory's range"""
        return (framenum >= self.beginFrame) and (framenum < self.endFrame)

    def getFrame(self, int framenum):
        """Returns absolute frame"""
        return self.pointData[framenum - self.beginFrame]

    def split(self, int framenum):
        cdef int split_idx = framenum - self.beginFrame
        t1 = Traj(self.beginFrame, self.name)
        t1.pointData = self.pointData[:split_idx]
        t2 = Traj(framenum, self.name)
        t2.pointData = self.pointData[split_idx:]
        return (t1, t2)

    def trimRight(self, int framenum):
        relFrame = self.getFrame(framenum)
        if relFrame <= 0:
            return None
        elif relFrame >= self.numFrames:
            return self
        else:
            t = Traj(self.beginFrame, self.name)
            t.pointData = self.pointData[:relFrame]
            return t

    def trimLeft(self, int framenum):
        relFrame = self.getFrame(framenum)
        if relFrame >= self.numFrames:
            return None
        elif relFrame <= 0:
            return self
        else:
            t = Traj(self.framenum, self.name)
            t.pointData = self.pointData[relFrame:]
            return t

    def predict(self, int idx=0):
        cdef int numpts = self.numFrames
        cdef int i
        if numpts == 0:     # No data, no predicton
            return None
        elif numpts == 1:   # One point, predict with itself
            return self.pointData[0]
        elif numpts == 2:   # Two points, fit a line
            result = []
            for i from 0 <= i < 3:
                x_2 = self.pointData[-2][i]
                x_1 = self.pointData[-1][i]
                result.append( x_2 + ( (x_1 - x_2) * (2 + idx) ) )
            return result
        else:               # Three or more points, fit a quadratic
            result = []
            for i from 0 <= i < 3:
                y1 = self.pointData[-3][i]
                y2 = self.pointData[-2][i]
                y3 = self.pointData[-1][i]
                a = y1
                b = -0.5*(y3 - 4*y2 + 3*y1)
                c = 0.5*(y3 - 2*y2 + y1)
                x = idx + 3
                result.append( c*x*x + b*x + a )
            return result

    def backPredict(self, int idx=0):
        cdef int numpts = self.endFrame - self.beginFrame
        if numpts == 0:     # No data, no predicton
            return None
        elif numpts == 1:   # One point, predict with itself
            return self.pointData[0]
        elif numpts == 2:   # Two points, fit a line
            fit = polyfit([1,2], array(self.pointData), 1)
            return list(polyval(fit, -idx))
        else:               # Three points, fit a quadratic
            fit = polyfit([1,2,3], array(self.pointData[-3:]), 2)
            return list(polyval(fit,-idx))

    def overlaps(self, other):
        "Returns true if there's a non empty intersection between frame ranges."
        return other.hasFrame(self.beginFrame) or other.hasFrame(self.endFrame)

    def contains(self, other):
        "Returns true if all frames of other are frames of self."



### Trajectorized data
######################

class TrajData(object):
    """Contains trajectories.

    Keeps a list of trajectories.
    Constructor can trajectorize from preprocessed RawData.

    """

    @property
    def numFrames(self):
        return \
            max([t.endFrame for t in self.trajs]) - \
            min([t.beginFrame for t in self.trajs])

    @property
    def maxFrame(self):
        return max([t.endFrame for t in self.trajs])

    @property
    def numTrajs(self):
        return len(self.trajs)

    def __init__(self):
        self.trajs = []
        self.filename = None
        self.framerate = 100.0
        self.changed = False
        self.trash = []

    def delete(self, traj):
        if traj in self.trajs:
            self.trash.append(traj)
            self.trajs.remove(traj)
        self.changed = True

    def undelete(self):
        if len(self.trash) > 0:
            self.trajs.append( self.trash.pop() )

    def splitTraj(self, traj, int framenum):
        if traj in self.trajs:
            t1,t2 = traj.split(framenum)
            self.trajs.remove(traj)
            self.trajs.extend([t1, t2])
        self.changed = True


    def cutRight(self, cutFrame):
        "Remove frames higher or equal than cutFrame"
        # Only keep trajectories that start before cutFrame
        self.trajs = [t for t in self.trajs if t.beginFrame < cutFrame ]
        # Cut trajectories that extend beyond numFram
        for t in self.trajs:
            if t.endFrame > cutFrame:
                newNumFrames = cutFrame - t.beginFrame
                t.pointData = t.pointData[:newNumFrames]
        self.changed = True

    def cutLeft(self, cutFrame):
        "Remove frames smaller than cutFrame"
        # Remove trajectories that end at cutFrame or before
        self.trajs = [t for t in self.trajs if t.endFrame >= cutFrame ]
        # Cut trajectories that start before cutFrame
        for t in self.trajs:
            if t.beginFrame < cutFrame:
                relFrame = cutFrame - t.beginFrame
                t.pointData = t.pointData[relFrame:]
                t.beginFrame = cutFrame
        # Make all trajectories start at cutFrame
        for t in self.trajs:
            t.beginFrame -= cutFrame
        self.changed = True
