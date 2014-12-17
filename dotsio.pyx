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

from trajdata import RawFrame, RawData, Traj, TrajData
import h5py
import numpy as np

# C3D
#####

def rawDataFromC3D(c3d, progress):
    cdef int i
    cdef int numFrames = c3d.header.numFrames()
    rd = RawData()
    rd.filename = c3d.filename
    rd.frameRate = c3d.header.frameRate

    for i from 0 <= i < numFrames:
        if i % 1000 == 0:
            progress.setValue( int(100.0*i / numFrames) )
            if progress.wasCanceled():
                return
        try:
            rd.frames.append( RawFrame(c3d.data[i,:,:]) )
            # if i == 20:
                # print (c3d.data[i,:,:])
        except:
            print('Error appending frame %d of %d' % (i, numFrames))
            break
    progress.setValue(100)
    return rd

# CSV
#####

def rawDataFromCSV(filename, progress):
    with open(filename) as pointsfile:
        framecount = 100
        rd = RawData()
        rd.filename = filename
        rd.frameRate = 100.0

        for line in pointsfile:
            line_array = line.split(',')
            if (line_array[0] == 'info'):
                if (line_array[1] == 'framecount'):
                    framecount = int(line_array[2])
            if (line_array[0] == 'frame'):
                i = len(rd.frames)
                if i % 1000 == 0:
                    progress.setValue( int(100.0*i / framecount) )
                if progress.wasCanceled():
                    return
                try:
                    rd.frames.append(readFrameFromArray(line_array, i))
                # if i == 20:
                    # print (c3d.data[i,:,:])
                except:
                   print('Error appending frame %d of %d' % (i, framecount))
                   break
        progress.setValue(100)
    
    if (not framecount):
        print ('Error frame count not found!')

    return rd
    
def readFrameFromArray(line, line_number):
    # This won't be used
    frame_id = int(line[1])
    frame_timestamp = float(line[2])
    rigidbody_count = int(line[3])
    counter = 4
    # read rigid bodies
    # for counter in range(rigidbody_count):
    counter = counter + rigidbody_count

    # This will
    frame_scale = 800
    marker_count = int(line[counter])
    if (marker_count == 0):
        print ('Warning, markers not found in the frame')
    data = np.zeros((marker_count, 4), dtype = np.float32)
    counter = counter + 1

    for i in range(marker_count):
        # 5 is the index for x in the first marker, and the other 5 is for the
        # number of variables (x, y, z, id, name) 
        #x
        data[i, 0] = np.float32(line[counter]) * frame_scale
        #y
        data[i, 2] = np.float32(line[counter + 1]) * frame_scale
        #z
        data[i, 1] = np.float32(line[counter + 2]) * -1 * frame_scale
        counter = counter + 5
    # if line_number == 0:
        # print (data)
    rf = RawFrame(data)
    return rf

# CSV 2
#######

def rawDataFromCSV2(filename, progress):
    with open(filename) as pointsfile:
        framecount = 100
        rd = RawData()
        rd.filename = filename
        rd.frameRate = 100.0
        counter = 0

        for line in pointsfile:
            counter = counter + 1
            # Line with the information
            if (counter == 1):
                line_array = line.split(',')
                framecount = int(line_array[11])
            # Line containing a frame
            if (counter > 7):
                line_array = line.split(',')
                i = len(rd.frames)
                if i % 1000 == 0:
                    progress.setValue( int(100.0*i / framecount) )
                if progress.wasCanceled():
                    return
                rd.frames.append(readFrameFromArray2(line_array, i))
        progress.setValue(100)
    
    if (not framecount):
        print ('Error frame count not found!')

    return rd
    
def readFrameFromArray2(line, line_number):
    # if (line_number == 0):
        # print (line)
    # This won't be used
    frame_id = int(line[0])
    frame_timestamp = float(line[1])

    # This will
    marker_count = 0
    points = []
    frame_scale = 800
    start = 2
    # read the indices in the array of the markers present in this frame
    for i in range(start, len(line)-1):
        if (line[i]):
            points.append(i)
    # if (line_number == 0):
        # print (points)
    data = np.zeros((len(points)/3, 4), dtype = np.float32)
    # divided by three because each marker has x, y, z
    for i in range(points[0], points[len(points)-1], 3):
        # if this cell is not empty
        if(line[i]):
            # print ("line: " + line[i])
            # print ("i: " + str(i))
            # print ("marker_count: " + str(marker_count))
            # read the next 3 elements
            #x
            data[marker_count, 0] = np.float32(line[i]) * frame_scale
            #y
            data[marker_count, 2] = np.float32(line[i + 1]) * frame_scale
            #z
            data[marker_count, 1] = np.float32(line[i + 2]) * -1 * frame_scale
            marker_count = marker_count + 1
    # if (line_number == 0):
        # print ("data: \n" + str(data))
    rf = RawFrame(data)
    return rf


# HDF5
######

def trajDataFromH5(filename, progress=None):
    """Read data from hdf5 file"""
    td = TrajData()
    td.filename = filename
    source = h5py.File(filename, 'r')
    group = source['trajectories']
    if not 'frame_rate' in group.attrs:
        print("Warning, no framerate specified. Setting to 100 FPS.")
        td.frameRate = 100.0
        td.changed = True
    else:
        td.frameRate = group.attrs['frame_rate']
    if not 'format_version' in group.attrs:
        print("Warning, no format version specified. Save to correct.")
        td.changed = True
    cdef int totalTrajs = len(list(group))        # Find number of trajectories
    for dsetName in list(group):
        dset = group[dsetName]
        tr = Traj(int(dset.attrs['begin_frame']), str(dset.attrs['name']))
        tr.pointData = dset.value.tolist()
        td.trajs.append(tr)
        if progress is not None:
            progress.setValue( int(100.0*td.numTrajs / totalTrajs) )
    source.close()
    del source
    if progress is not None:
        progress.setValue(100)
    return td

def trajToH5Dataset(h5group, idx, traj):
    shape = (traj.numFrames, 3)
    dset = h5group.create_dataset(
        'traj%d' % idx, shape, 'f', compression='gzip', data=traj.pointData)
    dset.attrs["name"] = traj.name
    dset.attrs["begin_frame"] = traj.beginFrame

def trajDataSaveH5(trajData, progress=None):
    f = h5py.File(trajData.filename, 'w')
    trajgroup = f.create_group('trajectories')
    trajgroup.attrs['format_version'] = 'dots 0'
    trajgroup.attrs['frame_rate'] = trajData.framerate
    cdef int idx = 0
    cdef int numTrajs = trajData.numTrajs
    for traj in trajData.trajs:
        idx += 1
        trajToH5Dataset(trajgroup, idx, traj)
        if progress is not None:
            progress.setValue( int(100.0 * idx / numTrajs) )
    f.flush()
    f.close()
    del f
    if progress is not None:
        progress.setValue(100)


# Raw to trajectorized data
###########################

def trajDataFromRawData(rdata, progress):
    cdef int numTrajs, frnum, min_k, numassig, numrem
    cdef float mindist, dx, dy, dz, sqDist
    td = TrajData()
    td.frameRate = rdata.frameRate

    # Trajectorization.
    # TODO: The following piece of code works, but I don't lke it.
    # A different class should trajectorize aRawData into a new TrajData.

    td.filename = rdata.filename + '.qtd'
    # Start with as many trajectories as points in the first frame
    currentTrajs = []
    trajIdx = 1
    numTrajs = len(rdata.frames[0].unusedIndices)
    for i in rdata.frames[0].unusedIndices:
        currentTrajs.append( Traj(0, "tr_%d" % trajIdx) )
        trajIdx += 1
        currentTrajs[-1].addPoint(rdata.frames[0].getSinglePoint(i))
    # For each frame
    frnum = -1
    for frame in rdata.frames[1:]:
        frnum += 1
        if frnum % 1000 == 0:
            progress.setValue( int(100.0*frnum / rdata.numFrames) )
            if progress.wasCanceled():
                return
        # Keep track of best point to add to each trajectory
        assignments = [-1 for t in currentTrajs]
        sqDists = [400 for t in currentTrajs]
        # All trajectory predictions
        predictions = [t.predict() for t in currentTrajs]
        # For each point in the frame
        for i in frame.unusedIndices:
            # Current point
            x0, y0, z0 = frame.getSinglePoint(i)
            # find closest trajectory
            min_dist = 2000000000.0
            min_k = -1
            for k in range(len(predictions)):
                # compute distance
                x,y,z = predictions[k]
                dx = x-x0
                dy = y-y0
                dz = z-z0
                sqDist = dx*dx + dy*dy + dz*dz
                # record distance if smaller
                if sqDist < min_dist:
                    min_dist = sqDist
                    min_k = k
            # assign
            if min_k != -1:
                if min_dist < sqDists[min_k]:
                    assignments[min_k] = i
                    sqDists[min_k] = min_dist
        # Apply assignments by adding to current trajectories.
        # Trajectories with no assignment are no longer current.
        numassig = 0
        numrem = 0
        no_longer_current = []
        for k in range(len(assignments)):
            if assignments[k] != -1:
                currentTrajs[k].addPoint(frame.getSinglePoint(assignments[k]))
                numassig += 1
                frame.unusedIndices.remove(assignments[k])
            else:
                no_longer_current.append(currentTrajs[k])
                numrem += 1
        # No longer current trajectories get removed from current and
        # added to self.trajs
        for tr in no_longer_current:
            td.trajs.append(tr)
            currentTrajs.remove(tr)
        # still unused points start new trajectories
        for i in frame.unusedIndices:
            newTraj = Traj(frnum, "tr_%d" % trajIdx)
            trajIdx += 1
            newTraj.addPoint(frame.getSinglePoint(i))
            currentTrajs.append(newTraj)
    # After processing all frames, still current trajectories get added to self.trajs
    td.trajs.extend(currentTrajs)
    progress.setValue(100)
    return td
