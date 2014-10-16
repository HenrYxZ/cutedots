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

from time import time
import cPickle as pickle
from trajdata import TrajData

# Body part coding
##################

# body part codes
bpname = [ 'UB', 'LB', 'Kn', 'Ft', 'El', 'Hn', 'tr', 'Hd' ]

# body part numbers
bpnum = dict([ (bpname[i],i) for i in range(len(bpname)) ])
numbp = len(bpname)  # number of body part types

# Human readable body parts
humanbp = { 'UB':'upper back', 'LB':'lower back', 'Kn':'knee',
            'Ft':'foot', 'El':'elbow', 'Hn':'hand', 'tr':'<none>', 'Hd':'Head' }

# Human readable body part sides
humanside = { 'L': 'left', 'R': 'right' }


cdef class ModelState:
    def __init__(ModelState self, data=None):
        self.data = data
        self.frame = 0
        self.traj = 0
        self.cam = CamState()
        self.anaglyph3d = 0
        self.bodyTrans = 0.7
        self.width = -1    # Reshape callback must set these
        self.height = -1
        self.aspect = -1.0

        if not self.data is None and self.data.numTrajs > 0:
            if not self.data.trajs[self.traj].hasFrame(self.frame):
                self.selNextTraj()
        self.trajMap = {}

    def currentTraj(ModelState self):
        if self.data is None or self.data.numTrajs == 0:
            return None
        if self.traj >= self.data.numTrajs:
            self.traj = self.data.numTrajs - 1
        if self.traj < 0:
            self.traj = 0
        return self.data.trajs[self.traj]

    # Transport
    def transBack(ModelState self, amount):
        self.frame = max(0, self.frame-amount)
    def transFwd(ModelState self, amount):
        self.frame = min(self.numFrames-1, self.frame+amount)
    def transStart(ModelState self):
        self.frame = 0
    def transEnd(ModelState self):
        self.frame = self.numFrames-1

    # Labeling
    def label(ModelState self, start, value):
        if self.data is None:
            return
        t = self.currentTraj()
        end = start + len(value)
        if t.name[start:end] == value:
            return
        t.name = t.name[:start] + value + t.name[end:4]
        self.data.changed = True
    def labelSubject(ModelState self, int subj):
        self.label(3, "%d" % subj)
    def labelSide(ModelState self, str side):
        self.label(2, side[0])
    def labelPartPrev(ModelState self):
        self.label(0, bpname[ max(0, self.bodyPartNum()-1) ])
    def labelPartNext(ModelState self):
        self.label(0, bpname[ min(numbp-1, self.bodyPartNum()+1) ])

    # Label parts
    def bodyPart(ModelState self):
        return self.currentTraj().name[:2]
    def side(ModelState self):
        return self.currentTraj().name[2]
    def subject(ModelState self):
        return self.currentTraj().name[3]

    # Info about labels
    def bodyPartNum(ModelState self):
        return bpnum[ self.bodyPart() ]
    def humanReadPart(ModelState self):
        n = self.currentTraj().name
        return humanside.get(n[2], '-') + ' ' + humanbp.get(n[:2], '-') + ' of ' + n[3]

    # Trajectory mapping
    def map(ModelState self, name, px, py, pz):
        """Maps trajectory names to positions.
        Also corrects for typical marker offsets just
        for visualization
        """

        x,y,z = px,py,pz
        part, side, subj = name[:2], name[2], name[3]
        # rise
        if part in ['Kn','UB']:
            z += 50
        # move back
        if part in ['Kn','Ft']:
            x = x-50 if subj=='1' else x+50
        # move forward
        if part in ['UB']:
            x = x+100 if subj=='1' else x-100
        # LB
        if part == 'LB':
            # move inside
            if subj == '1':
                y = y-80 if side == 'L' else y+80
            else:
                y = y+80 if side == 'L' else y-80
            # lower
            z -= 100
            # move forward
            x = x+20 if subj=='1' else x-20
        # UB
        if part == 'UB':
            # Move outside
            if subj == '1':
                y = y+70 if side == 'L' else y-70
            else:
                y = y-70 if side == 'L' else y+70

        self.trajMap[name] = (x, y, z)

    def selNextTraj(ModelState self):
        if self.data is None:
            return
        cdef int start = self.traj
        self.traj = (self.traj+1) % self.data.numTrajs
        while self.traj != start:
            if self.currentTraj().hasFrame(self.frame):
                break
            self.traj = (self.traj+1) % self.data.numTrajs

    def selPrevTraj(ModelState self):
        if self.data is None:
            return
        cdef int start = self.traj
        self.traj = (self.data.numTrajs + self.traj - 1) % self.data.numTrajs
        while self.traj != start:
            if self.currentTraj().hasFrame(self.frame):
                break
            self.traj = (self.data.numTrajs+self.traj-1) % self.data.numTrajs

    def currentIsUnlabeled(ModelState self):
        if self.bodyPart() != 'Hd':
            return  not self.bodyPart() in bpname[:6] \
                or not self.side() in ['L', 'R'] \
                or not self.subject() in ['1', '2']
        else:
            return not self.subject() in ['1', '2']

    def selNextUnlabTraj(ModelState self):
        if self.data is None:
            return
        cdef int start = self.traj
        self.traj = (self.traj+1) % self.data.numTrajs
        while self.traj != start:
            if self.currentIsUnlabeled():
                self.frame = self.currentTraj().beginFrame
                break
            self.traj = (self.traj+1) % self.data.numTrajs

    def selPrevUnlabTraj(ModelState self):
        if self.data is None:
            return
        cdef int start = self.traj
        self.traj = (self.data.numTrajs+self.traj-1) % self.data.numTrajs
        while self.traj != start:
            if self.currentIsUnlabeled():
                self.frame = self.currentTraj().beginFrame
                break
            self.traj = (self.data.numTrajs+self.traj-1) % self.data.numTrajs

    def selNextBreak(ModelState self):
        if self.data is None:
            return
        currentTrajs = [x for x in self.data.trajs if x.hasFrame(self.frame)]
        if len(currentTrajs) > 0:
            byEndFrame = sorted(currentTrajs, key=lambda x: x.endFrame)
            self.frame = byEndFrame[0].endFrame

    def deleteCurrent(ModelState self):
        if self.data is None:
            return
        self.data.delete( self.currentTraj() )
        self.traj = self.traj % self.data.numTrajs

    def undelete(ModelState self):
        if self.data is None:
            return
        self.data.undelete()

    def deleteUnassigned(self):
        if self.data is None:
            return
        deleteList = []
        for traj in self.data.trajs:
            bp, side, subj = traj.name[:2], traj.name[2], traj.name[3]
            if (bp == 'tr') or (not (bp in bpname)) or \
                (not (subj in ['1', '2'])) or \
                ((bp != 'Hd') and (not (side in ['L', 'R']))):
                    deleteList.append(traj)
        if len(deleteList) > 0:
            for traj in deleteList:
                self.data.trajs.remove(traj)
            self.data.changed = True

    def cutRight(self):
        if self.frame > 0:
            self.data.cutRight(self.frame)
            self.frame = self.data.numFrames - 1
        else:
            self.data = None
            self.traj = 0

    def cutLeft(self):
        if self.frame < self.data.numFrames:
            self.data.cutLeft(self.frame)
            self.frame = 0
        else:
            self.data = None
            self.traj = 0
