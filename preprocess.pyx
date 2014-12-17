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
import dotsio

def ppC3D(filename, progress):
    if progress.wasCanceled():
        return
    # Read C3D data
    progress.setLabelText( "Reading file..." )
    c3d = c3dformat.C3d(filename)
    if progress.wasCanceled():
        return
    # Read raw data
    progress.setLabelText( "Extracting non-trajectorized data..." )
    rd = dotsio.rawDataFromC3D(c3d, progress)
    if progress.wasCanceled():
        return
    # Remove duplicates
    progress.setLabelText( "Merging close points..." )
    rd.joinClosePoints(progress)
    if progress.wasCanceled():
        return
    # Initial trajectorization
    progress.setLabelText( "Initial trajectorization..." )
    td = dotsio.trajDataFromRawData(rd, progress)
    if progress.wasCanceled():
        return
    # Writing
    progress.setLabelText( "Saving trajectorized data..." )
    dotsio.trajDataSaveH5(td, progress)
    return td.filename

def ppCSV(filename, progress):
    if progress.wasCanceled():
        return
    # Read raw data
    progress.setLabelText( "Extracting non-trajectorized data..." )
    rd = dotsio.rawDataFromCSV(filename, progress)
    if progress.wasCanceled():
        return
    # Remove duplicates
    progress.setLabelText( "Merging close points..." )
    rd.joinClosePoints(progress)
    if progress.wasCanceled():
        return
    # Initial trajectorization
    progress.setLabelText( "Initial trajectorization..." )
    td = dotsio.trajDataFromRawData(rd, progress)
    if progress.wasCanceled():
        return
    # Writing
    progress.setLabelText( "Saving trajectorized data..." )
    dotsio.trajDataSaveH5(td, progress)
    return td.filename

def ppCSV2(filename, progress):
    if progress.wasCanceled():
        return
    # Read raw data
    progress.setLabelText( "Extracting non-trajectorized data..." )
    rd = dotsio.rawDataFromCSV2(filename, progress)
    if progress.wasCanceled():
        return
    # Remove duplicates
    progress.setLabelText( "Merging close points..." )
    rd.joinClosePoints(progress)
    if progress.wasCanceled():
        return
    # Initial trajectorization
    progress.setLabelText( "Initial trajectorization..." )
    td = dotsio.trajDataFromRawData(rd, progress)
    if progress.wasCanceled():
        return
    # Writing
    progress.setLabelText( "Saving trajectorized data..." )
    dotsio.trajDataSaveH5(td, progress)
    return td.filename