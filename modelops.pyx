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

from trajdata import Traj
import numpy as np

class ProcessingException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return str(self.msg)


def guessSideAndSubject(data):
    for t in data.trajs:
        n = t.name
        if n[:3] == 'tr_':
            avg = t.average()
            if avg[0] < 0.0:
                if avg[1] < 0.0: t.name = 'trR1'
                else:            t.name = 'trL1'
            else:
                if avg[1] < 0.0: t.name = 'trL2'
                else:            t.name = 'trR2'
    data.changed = True


def sortTrajs(data):
    data.trajs = sorted(data.trajs, key=lambda t: t.pointData[0][0])

def sortTrajsSlow(data):
    data.trajs = sorted(data.trajs, key=lambda t: t.averageX())

def rotate90Deg(data, progress):
    cdef int trajNum = 0
    cdef int numTrajs = data.numTrajs
    for traj in data.trajs:
        for i in range(traj.numFrames):
            x, y, z = traj.pointData[i]
            traj.pointData[i] = -y, x, z
        trajNum += 1
        progress.setValue(int(100.0 * trajNum / numTrajs))
    data.changed = True
    progress.setValue(100)


def swapSubjects(data, progress):
    for traj in data.trajs:
        n = traj.name
        if n[3] == '1':
            traj.name = n[:3] + '2'
        elif n[3] == '2':
            traj.name = n[:3] + '1'
    rotate90Deg(data, progress)
    rotate90Deg(data, progress)
 

def averageSameNameTrajectories(data, progress):
    cdef float sx, sy, sz
    cdef int count, i
    # index trajectories by name
    indTrajs = {}
    heads = [] 
    for traj in data.trajs:
        if traj.name[:2] == 'Hd':
            heads.append(traj)
        else:
            indTrajs.setdefault(traj.name, []).append(traj)
    newTrajs = [] # List of new trajedies
    progress.setValue(0)
    # average each list of trajectories
    nameCount = 0
    for name, trajs in indTrajs.items():
        # find being and end frame
        beginFrame = min([t.beginFrame for t in trajs])
        endFrame = max([t.endFrame for t in trajs])
        current = None
        for i from beginFrame <= i < endFrame:
            # average points in frame
            sx = sy = sz = 0.0
            count = 0
            for x, y, z in  [t.getFrame(i) for t in trajs if t.hasFrame(i)]:
                sx += x
                sy += y
                sz += z
                count += 1
            if count >= 1:  # Append point average
                if current is None:
                    current = Traj(i, name)
                    newTrajs.append(current)
                current.pointData.append([sx/count, sy/count, sz/count])
            else:           # Close current trajectory
                current = None
        nameCount += 1
        progress.setValue(int(100.0*nameCount / len(indTrajs.keys())))
    newTrajs.extend(heads)
    data.trajs = newTrajs
    progress.setValue(100)
    data.changed = True

def fillGaps(data, maxGapTime, maxSampleTime, progress):
    """Joins same name trajectories together by interpolating gaps between them.

    Arguments:
    data          -- TrajData object from which to take trajectories.
    maxGapTime    -- Two trajectories will only be connected if the gap to interpolate
                     between them is this amount of time or less.
    maxSampleTime -- Consider at most this time from the end of a trajectory and at most
                     the same time from the begining of the next, for quadratic
                     interpolation at gap points.
    """

    maxGap = int(data.framerate * maxGapTime)
    maxSample = int(data.framerate * maxSampleTime)
    # index trajectories by name
    indTrajs = {}
    heads = [] 
    for traj in data.trajs:
        if traj.name[:2] == 'Hd':
            heads.append(traj)
        else:
            indTrajs.setdefault(traj.name, []).append(traj)
    deleteList = []
    for name, trajs in indTrajs.items():
        # sort trajs by first frame
        trajs = sorted(trajs, key=lambda t: t.beginFrame)
        t0 = trajs[0]
        for t in trajs[1:]:
            gap = t.beginFrame - t0.endFrame    # find gap length
            if gap < 0:
                raise ProcessingException("Same name trajectories overlap. " + \
                    "Make sure trajectories are labeled correctly and use " + \
                    "'Average same name trajectories' operation to correct minor overlaps")
            elif gap == 0:
                t0.pointData.extend(t.pointData)
                deleteList.append(t)
            elif gap > 0 and gap <= maxGap:
                # Last points of first trajectory
                fromPts = t0.pointData[-maxSample:]
                # First points of second trajectory
                toPts = t.pointData[:maxSample]
                # Data for model fit
                trainData = fromPts + toPts
                fromLen, toLen = len(fromPts), len(toPts)
                fromTime = range(fromLen)
                gapTime = range(fromLen, fromLen + gap)
                toTime = range(fromLen + gap, fromLen + gap + toLen)
                time = fromTime + toTime
                x = [point[0] for point in trainData]
                y = [point[1] for point in trainData]
                z = [point[2] for point in trainData]
                order = 2
                if len(time) <=2:
                    order = 1
                fitx = np.polyfit(time, x, 2)
                fity = np.polyfit(time, y, 2)
                fitz = np.polyfit(time, z, 2)
                gapx = np.polyval(fitx, gapTime)
                gapy = np.polyval(fity, gapTime)
                gapz = np.polyval(fitz, gapTime)
                gapPoints = [ [gapx[fr], gapy[fr], gapz[fr]]
                               for fr in range(len(gapTime)) ]
                t0.pointData.extend(gapPoints)
                t0.pointData.extend(t.pointData)
                deleteList.append(t)
            else:
                t0 = t
    for t in deleteList:
        data.trajs.remove(t)
    data.changed = True
