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

import pylab as pl
import numpy as np
import scipy as sp


def plotContinuity(data):
    # Organize curves by label
    trajectories = dict()
    for traj in data.trajs:
        trajectories.setdefault(traj.name, []).append(traj)
    # Plot
    y = 0
    ticksy = []
    ticksn = []
    for name, trajs in trajectories.items():
        y -= 10
        ticksy.append(y)
        ticksn.append(name)
        for t in trajs:
            t0, t1 = t.beginFrame / data.framerate, (t.endFrame-1) / data.framerate
            pl.plot([t0, t1], [y, y])
    pl.yticks(ticksy, ticksn)
    pl.show()

def energy(trajectories):
    "Returns aggregated energy signal for all markers"

    maxFrame = max([t.endFrame for t in trajectories])
    energy = np.zeros(maxFrame-1) # Speed results
    for t in trajectories:
        if t.numFrames >= 2:
            trajSpeed = np.subtract(100*np.array(t.pointData[1:]),
                                    100*np.array(t.pointData[:-1]))  # find speed
            trajEnergy = np.sum(trajSpeed**2,1)
            energy[t.beginFrame:t.endFrame-1] += trajEnergy     # store energy
    return energy


def trajsBySubj(trajectories):
    # classify trajectories by subject
    trajs = {}
    for t in trajectories:
        trajs.setdefault(t.name[3], []).append(t)
    return trajs

def plotSubjectFunc(data, func, subplot=True):
    "Makes a plot of a time function by subject"

    trajs = trajsBySubj(data.trajs)
    numSubjs = len(trajs)
    subjNum = 1

    for subject, trajList in trajs.items():
        f = func(trajList)
        t = np.arange(len(f)) / data.framerate
        if subplot:
            pl.subplot(numSubjs, 1, subjNum)
            pl.ylabel("Subject %s" % subject)
        subjNum += 1
        pl.plot(t, f)
 

def plotEnergy(data):
    plotSubjectFunc(data, energy)
    pl.show()

def ortoSpeeds(trajectories):
    "Returns aggregation of speeds for each axis."

    maxFrame = max([t.endFrame for t in trajectories])
    speed = np.zeros((maxFrame-1,3)) # Speed results
    for t in trajectories:
        if t.numFrames >= 2:
            trajSpeed = np.subtract(t.pointData[1:], t.pointData[:-1])  # find speed
            speed[t.beginFrame:t.endFrame-1] += trajSpeed  # add traj. speed to global
    return speed

def plotSpeed(data):
    # classify trajectories by subject
    trajs = trajsBySubj(data.trajs)

    axes = [ ('horizontal', 0), ('sides', 1), ('vertical', 2) ]
    plotNum = 1
    subjSpeeds = [ (subj, ortoSpeeds(trajs)) for (subj,trajs) in trajs.items() ]
    numPlots = len(axes)

    for title, axis in axes:
        pl.subplot(numPlots,1,plotNum)
        plotNum += 1
        for subj, speed in subjSpeeds:
            t = np.arange(len(speed)) / data.framerate
            pl.plot(t, speed[:,axis])
            pl.ylabel(title)
    pl.show()


def projectedSpeed(trajectories, direction):
    "Returns aggregation of projected speeds onto the given (unitary) direction"

    maxFrame = max([t.endFrame for t in trajectories])
    speed = np.zeros(maxFrame-1) # Speed results
    for t in trajectories:
        if t.numFrames >= 2:
            trajSpeed = np.subtract(t.pointData[1:], t.pointData[:-1])  # find speed
            projSpeed = np.sum(trajSpeed * direction, 1)
            speed[t.beginFrame:t.endFrame-1] += projSpeed  # store energy
    return speed

def plotLengthHistogram(data):
    pl.hist([t.numFrames / data.framerate for t in data.trajs])
    pl.xlabel('length (secs)')

def plotLengthHistAll(data):
    plotLengthHistogram(data)
    pl.show()
