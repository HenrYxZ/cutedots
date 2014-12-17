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


# Display size
DISPLAY_WIDTH = 1024
DISPLAY_HEIGHT = 600

from PySide import QtGui, QtOpenGL, QtCore
from PySide.QtCore import Qt
import dotsdisplay as dots
import trajdata
import modelops
import dotsio
import analysis
import preprocess
import os
import sys
import time
import traceback

about_text = """
Copyright 2012 Esteban Hurtado

Cutedots is distributed under the terms of the Reciprocal Public License 1.5.

Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.

You should have received a copy of the Reciprocal Public License along with Cutedots. If not, see <http://opensource.org/licenses/rpl-1.5>.
"""

keybindings = [
    ('Navigate frames', [
            ('next/prev. frame', 'arrows'),
            ]),
    ('View', [
            ('rotate','a d'),
            ('zoom', 's w'),
            ('reset view', '0'),
            ]),
    ('Curve selection', [
            ('select prev, next', 'f g'),
            ('unlabeled prev, next', 'r t'),
            ('select next break', 'b'),
            ]),
    ('Labeling', [
            ('subject', '1 2'),
            ('left, right', 'n m'),
            ('part prev, next', ', .'),
            ]),
    ('Modify', [
            ('split trajectory (cut)', 'c'),
            ('delete trajectory', 'z'),
            ('undelete trajectory', 'u'),
#            ('delete unassigned', 'k'),
            ]),
    ('Other', [
            ('super cool anaglyph :)', '3'),
            ])
]


def helptext():
    text = "<table>\n"
    for title, bindings in keybindings:
        text += '  <tr><th colspan="4" align="left">%s</th></tr>\n' % title
        for b in bindings:
            text += ('  <tr><td width="20px"></td><th align="left">%s</th>' + \
                         '<td width="20px"></td><td>%s</td></tr>\n') % b
        text += '<tr><td colspan=4></td></tr>\n'
    text += '</table>\n'
    return text


class GLWidget(QtOpenGL.QGLWidget):
    def __init__(self, parent):
        QtOpenGL.QGLWidget.__init__(self, parent)
        self.parentWindow = parent
        self.setMinimumSize(640, 480)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.updateGL)
        self.timer.timeout.connect(self.parent().updateFrame)
        self.timer.timeout.connect(self.parent().updateStatus)
        self.setFocusPolicy(Qt.StrongFocus)
        self.dots = dots.DotsDisplay(self.getTextRenderFunc())
        

    def getTextRenderFunc(self):
        return lambda s: self.renderText(0.0, 0.0, 0.0, s)

    def paintGL(self)        : self.dots.display()
    def resizeGL(self, w, h) : self.dots.reshape(w, h)
    def initializeGL(self)   : self.dots.setupViewport(DISPLAY_WIDTH, DISPLAY_HEIGHT)

    def togglePlay(self):
        # TODO
        return
        # dots.disp_data().toggle_play()
        # if dots.disp_data().playing == 1:
        #     self.timer.start(25)
        # else:
        #     self.timer.stop()
 
    def keyPressEvent(self, e):
        key = e.key()
        model = self.dots.model
        # Change frame
        if key == Qt.Key_Right:     model.frame = min(model.frame+1, model.data.numFrames-1)
        elif key == Qt.Key_Left:  model.frame = max(model.frame-1, 0)
        # Rotation
        elif   key == Qt.Key_A:    model.cam.rotation -= 2
        elif key == Qt.Key_D:      model.cam.rotation += 2
        # Zoom
        elif key == Qt.Key_S:      model.cam.py -= 25  
        elif key == Qt.Key_W:      model.cam.py += 25
        # Curve selection
        elif key == Qt.Key_R:      model.selPrevUnlabTraj()
        elif key == Qt.Key_T:      model.selNextUnlabTraj()
        elif key == Qt.Key_F:      model.selPrevTraj()
        elif key == Qt.Key_G:      model.selNextTraj()
        elif key == Qt.Key_B:      model.selNextBreak()
        # Anaglyph
        elif key == Qt.Key_3:      model.anaglyph3d=not model.anaglyph3d
        # Cam reset
        elif key == Qt.Key_0:      model.cam.reset()
        # Curve labeling
        elif key == Qt.Key_1:      model.labelSubject(1)
        elif key == Qt.Key_2:      model.labelSubject(2)
        elif key == Qt.Key_N:      model.labelSide('L')
        elif key == Qt.Key_M:      model.labelSide('R')
        elif key == Qt.Key_Comma:  model.labelPartPrev()
        elif key == Qt.Key_Period: model.labelPartNext()
        # Curve splitting (cut)
        elif key == Qt.Key_C:      model.data.splitTraj(model.currentTraj(), model.frame) 
       # Curve deletion
        elif key == Qt.Key_Z:      model.deleteCurrent()
        elif key == Qt.Key_U:      model.undelete()
#        elif key == Qt.Key_K:      model.deleteUnassigned()
        self.updateGL()
        self.parentWindow.updateStatus()
        self.parentWindow.updateFrame()

class FrameSpin(QtGui.QSpinBox):
    def __init__(self, parent):
        QtGui.QSpinBox.__init__(self, parent)
        self.setWrapping(True)
        self.setAccelerated(True)

    def textFromValue(self, val):
        mins = val / 6000
        secs = (val / 100) % 60
        frames = val % 100
        return "%02d':%02d\".%02d" % (mins, secs, frames)


class TransportBar(QtGui.QWidget):
    def __init__(self, parent, glwidget):
        QtGui.QWidget.__init__(self, parent)
        self.mainWindow = parent
        self.gl = glwidget
        self.playing = False
        self.refFrame = 0
        self.refTime = 0.0
        self.playSpeed = 1.0
        self.initUI()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update)
        self.playButton.released.connect(self.togglePlay)

    def initUI(self):
        # Get play/pause icons
        standardIcon = QtGui.qApp.style().standardIcon
        self.iconPlay = standardIcon(QtGui.QStyle.SP_MediaPlay)
        self.iconPause = standardIcon(QtGui.QStyle.SP_MediaPause)
        # Widgets
        self.playButton = QtGui.QPushButton(self.iconPlay, '', self)
        s1 = self.speed1xButton = QtGui.QToolButton(self)
        s1.setText('1x')
        s2 = self.speed4xButton = QtGui.QToolButton(self)
        s2.setText('4x')
        s3 = self.speed16xButton = QtGui.QToolButton(self)
        s3.setText('16x')
        self.frameSpin = FrameSpin(self)
        self.frameSlider = QtGui.QSlider(Qt.Horizontal, self)
        self.labelTotal = QtGui.QLabel(self)
        # Button group
        buttonGroup = QtGui.QButtonGroup(self)
        for button in [s1, s2, s3]:
            button.setCheckable(True)
            buttonGroup.addButton(button)
        s1.setChecked(True)
        # Layouts
        layout = QtGui.QHBoxLayout(self)
        layout.addWidget(self.playButton)
        layout.addWidget(self.speed1xButton)
        layout.addWidget(self.speed4xButton)
        layout.addWidget(self.speed16xButton)
        layout.addWidget(self.frameSpin)
        layout.addWidget(self.frameSlider)
        layout.addWidget(self.labelTotal)
        # Signals
        self.frameSpin.valueChanged.connect(self.frameSlider.setValue)
        self.frameSpin.valueChanged.connect(self.parent().setFrame)
        self.frameSlider.valueChanged.connect(self.frameSpin.setValue)
        self.speed1xButton.clicked.connect(self.setSpeed)
        self.speed4xButton.clicked.connect(self.setSpeed)
        self.speed16xButton.clicked.connect(self.setSpeed)

    def setReference(method):
        def deco(self, *args, **kargs):
            self.refFrame = self.gl.dots.model.frame
            self.refTime = time.time()
            method(self, *args, **kargs)
        return deco

    @setReference
    def play(self):
        self.playing = True
        self.playButton.setIcon(self.iconPause)
        self.timer.start(25)

    @setReference
    def setSpeed(self):
        if self.speed1xButton.isChecked():
            self.playSpeed = 1.0
        elif self.speed4xButton.isChecked():
            self.playSpeed = 4.0
        elif self.speed16xButton.isChecked():
            self.playSpeed = 16.0

    def stop(self):
        self.timer.stop()
        self.playing = False
        self.playButton.setIcon(self.iconPlay)

    def togglePlay(self):
        if self.playing: self.stop()
        else:            self.play()

    def update(self):
        t = time.time()
        model = self.gl.dots.model
        fr = model.data.frameRate * self.playSpeed
        model.frame = int(self.refFrame + (t - self.refTime) * fr) % model.data.numFrames
        self.gl.updateGL()
        self.mainWindow.updateFrame()
        self.mainWindow.updateStatus()

    def setNumFrames(self, numFrames):
        self.frameSpin.setRange(0, numFrames-1)
        self.frameSlider.setRange(0, numFrames-1)
        self.labelTotal.setText(self.frameSpin.textFromValue(numFrames-1))

    def setFrame(self, numFrame):
        spinBlocked = self.frameSpin.blockSignals(True)
        sliderBlocked = self.frameSlider.blockSignals(True)
        self.frameSpin.setValue(numFrame)
        self.frameSlider.setValue(numFrame)
        self.frameSpin.blockSignals(spinBlocked)
        self.frameSlider.blockSignals(sliderBlocked)

class CuteDotsMainWindow(QtGui.QMainWindow):

    @property
    def model(self):
        return self.gl.dots.model

    @property
    def data(self):
        return self.gl.dots.model.data

    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.initUI()

    def initUI(self):
        self.initMenu()
        # Widgets
        self.mainWidget = QtGui.QWidget(self)        # Main widget
        self.gl = GLWidget(self)                     # 3D view
        self.transport = TransportBar(self, self.gl) # TransportBar
        self.help = QtGui.QLabel(helptext(), self)   # Keyboard help text
        # Layout
        self.setCentralWidget(self.mainWidget)
        self.mainLayout = QtGui.QVBoxLayout(self.mainWidget)
        self.hwidget = QtGui.QWidget(self)
        self.hlayout = QtGui.QHBoxLayout(self.hwidget)
        self.hlayout.addWidget(self.gl, 1)
        self.hlayout.addWidget(self.help)
        self.mainLayout.addWidget(self.hwidget, 1)
        self.mainLayout.addWidget(self.transport)
        self.setGeometry(300, 300, 300, 300)
        self.setWindowTitle("Cute dots")
        self.updateStatus()

    # Decorators
    ############

    def ifDataLoaded(method):
        "Decorator that executes method only if data is loaded"
        def deco(self, *args, **kargs):
            if self.data is None:
                return
            method(self, *args, **kargs)
        return deco

    def warnIfNoDataLoaded(method):
        "Decorator that throws a warning if no data is already loaded."
        def deco(self, *args, **kargs):
            if self.data is None:
                QtGui.QMessageBox.warning(self, 'Error', 'No data loaded yet.')
                return
            method(self, *args, **kargs)
        return deco

    def updateDisplay(method):
        "Decorator that updates GL display after calling method."
        def deco(self, *args, **kargs):
            method(self, *args, **kargs)
            self.gl.updateGL()
        return deco

    def updateStatus(method):
        "Updates status after calling method."
        def deco(self, *args, **kargs):
            method(self, *args, **kargs)
            self.updateStatus()
        return deco

    def updateNumFrames(method):
        def deco(self, *args, **kargs):
            method(self, *args, **kargs)
            self.transport.setNumFrames(self.data.numFrames)
        return deco
        

    # End decorators #

    def mkAction(self, text, tip, handler, shortCut=None):
        act = QtGui.QAction(text, self)
        if not shortCut is None:
            act.setShortcut(shortCut)
        act.setStatusTip(tip)
        act.triggered.connect(handler)
        return act

    def mkMenu(self, menubar, name, actions):
        menu = menubar.addMenu(name)
        for action in actions:
            if action is None:
                menu.addSeparator()
            else:
                menu.addAction(action)

    def mkProgress(self, task_descript, total=100, cancelbtn=""):
        progress = QtGui.QProgressDialog(
            task_descript, cancelbtn, 0, total, self)
        progress.setWindowModality(Qt.WindowModal)
        progress.setWindowTitle("Cutedots")
        if cancelbtn == "":
            progress.setCancelButton(None)
        return progress

    def initMenu(self):
        self.statusBar()
        # File
        loadAction = self.mkAction('&Open motion data file...',
            'Open pickled trajectory data or import C3D file',
            self.openDataFile, 'Ctrl+O')
        importC3DAction = self.mkAction('&Import C3D file...',
            'Read and trajectorize data from C3D file (creates new CuteDots file on disk)',
            self.importC3D, 'Ctrl+I')
        importCSVAction = self.mkAction('&Import CSV file...',
            'Read and trajectorize data from CSV file (creates new CuteDots file on disk)',
            self.importCSV)
        importCSV2Action = self.mkAction('&Import CSV2 file...',
            'Read from CSV file using Motive <1.7', self.importCSV2)
        saveAction = self.mkAction('&Save', 'Save to HDF5 file (.qtd).',
             self.saveFile, 'Ctrl+S')
        saveAsAction = self.mkAction('Save &as...', 'Save to a different file',
             self.saveFileAs)
        saveSeqAction = self.mkAction('Save image &sequence',
            'Save frame to frame image sequence of realtime motion' + \
            'from current view', self.saveSequence)
        exitAction = self.mkAction(
            'Save and e&xit', 'Save and exit program', self.close, 'Ctrl+Q')
        # Operations
        guessAction = self.mkAction('Guess side and subject',
            'Label side and subject based on horizontal quadrant',
            self.guessSideAndSubject)
        rotateAction = self.mkAction('Rotate 90 degrees',
            'Rotate points 90 degrees in horizontal plane', self.rotate)
        swapAction = self.mkAction('Swap subjects',
            'Exchange subjects 1 and 2', self.swap)
        averageTrajsAction = self.mkAction('Average same name trajectories',
            'Makes sure only one point exists at each frame for each trajectory name',
            self.averageTrajs)
        fillSmallGapsAction = self.mkAction('Fill small gaps in trajectories',
            'Interpolates up to 0.2 seconds of missing frames in a trajectory.',
             self.fillSmallGaps)
        fillAllGapsAction = self.mkAction('Fill all gaps in trajectories',
            'Interpolates all missing frames in trajectories.', self.fillAllGaps)
        cutRightAction = self.mkAction('Cut from this frame', '', self.cutRight)
        cutLeftAction = self.mkAction('Cut before this frame', '', self.cutLeft)
        removeUnassignedAction = self.mkAction('Remove unassigned trajectories.',
            'Removes trajectories that are not completeley assigned.', self.removeUnassigned)
        sortAction = self.mkAction('Sort trajectories by X coordinate',
            'Improves traject. traversal order for labeling', self.slowSort)
        # Analysis
        plotContAction = self.mkAction('Plot continuity',
            'Plot frame span of each trajectory', self.plotContinuity)
        plotEnergyAction = self.mkAction('Plot energy',
            'Make global energy plot', self.plotEnergy)
        plotSpeedAction = self.mkAction('Plot speed', 'Plots speed for each axis.',
            self.plotSpeed)
        lengthHistAction = self.mkAction('Plot histogram of trajectory lengths', '',
            self.plotLengthHistogram)
        # Other
        exitNoSaveAction = self.mkAction(
            'Exit without saving', 'Exit the program', self.exitNoSave)
        # Help
        aboutAction = self.mkAction('&About', 'About', self.displayAboutDlg)
        # Menu bar
        menubar = self.menuBar()
        self.mkMenu(menubar, '&File', [loadAction, importC3DAction, importCSVAction, importCSV2Action,
            None, saveAction, saveAsAction, saveSeqAction, None, exitAction, exitNoSaveAction])
        self.mkMenu(menubar, '&Operations', [rotateAction, swapAction, guessAction, None,
            averageTrajsAction, fillSmallGapsAction, fillAllGapsAction, None,
            cutRightAction, cutLeftAction, None, removeUnassignedAction, None, sortAction])
        self.mkMenu(menubar, '&Analysis', [plotContAction, plotEnergyAction, plotSpeedAction,
            lengthHistAction])
        self.mkMenu(menubar, '&Help', [aboutAction])

    # Actions
    #########

    def displayAboutDlg(self):
        QtGui.QMessageBox.about(self, 'About', about_text)

    def openDataFile(self):
        fn,x = QtGui.QFileDialog.getOpenFileName(None, 'Open a CuteDots motion data file')
        if fn == '':
            return
        self.loadDataFile(fn)

    @updateStatus
    @updateNumFrames
    def loadDataFile(self, fn):
        if self.data is not None and self.data.changed:
            self.saveFile()
        progress = self.mkProgress("Loading C3D...", 100, "Cancel")
        self.gl.dots.loadData(fn, progress)
        progress.close()
        modelops.sortTrajs(self.data)
        del progress
        self.setWindowTitle("Cute dots - " + fn)

    def importC3D(self):
        fn,x = QtGui.QFileDialog.getOpenFileName(None, 'Import a C3D data file')
        if fn == '':
            return
        progress = self.mkProgress("Importing...")
        qtdFn = preprocess.ppC3D(fn, progress)
        progress.close()
        del progress
        self.loadDataFile(qtdFn)

    def importCSV(self):
        fn,x = QtGui.QFileDialog.getOpenFileName(None, 'Import a CSV data file')
        if fn == '':
            return
        progress = self.mkProgress("Importing...")
        qtdFn = preprocess.ppCSV(fn, progress)
        progress.close()
        del progress
        self.loadDataFile(qtdFn)

    def importCSV2(self):
        fn,x = QtGui.QFileDialog.getOpenFileName(None, 'Import a CSV2 data file')
        if fn == '':
            return
        progress = self.mkProgress("Importing...")
        qtdFn = preprocess.ppCSV2(fn, progress)
        progress.close()
        del progress
        self.loadDataFile(qtdFn)

    @ifDataLoaded
    def saveFile(self):
        if self.data.changed:
            progress = self.mkProgress("Saving...", 100, "Cancel")
            dotsio.trajDataSaveH5(self.data, progress)
            self.data.changed = False
            del progress

    @ifDataLoaded
    def saveFileAs(self):
        fn,x = QtGui.QFileDialog.getOpenFileName(None, 'Save motion data file')
        if fn == '':
            return
        self.data.changed = True
        self.data.filename = fn 
        self.saveFile()


    @warnIfNoDataLoaded
    def saveSequence(self):
        folder = QtGui.QFileDialog.getExistingDirectory(
            self, 'Choose a folder to save image files for the video sequence')
        if folder == '':
            return
        progress = self.mkProgress(
            "Generating images...", self.data.numFrames, "Stop")
        period = 100. / 30.
        step = 0
        frameNum = 0
        while frameNum < 6000:
            self.model.frame = frameNum
            self.gl.updateGL()
            self.updateStatus()
            self.gl.grabFrameBuffer().save(
                os.path.join(folder,'frame%06d.jpg' % step), None, 90)
            step += 1
            frameNum = int(step * period)
            if progress.wasCanceled():
                break
            progress.setValue(frameNum)
        progress.setValue(self.model.numFrames)

    @warnIfNoDataLoaded
    @updateDisplay
    def guessSideAndSubject(self):
        modelops.guessSideAndSubject(self.data)

    @warnIfNoDataLoaded
    @updateDisplay
    def rotate(self):
        progress = self.mkProgress("Rotating 90 degrees...")
        modelops.rotate90Deg(self.data, progress)
        progress.close()
        del progress

    @warnIfNoDataLoaded
    @updateDisplay
    def swap(self):
        progress = self.mkProgress("Swapping subjects...")
        modelops.swapSubjects(self.data, progress)

    @warnIfNoDataLoaded
    @updateDisplay
    @updateStatus
    def averageTrajs(self):
        progress = self.mkProgress("Averaging same name trajectories...")
        modelops.averageSameNameTrajectories(self.data, progress)

    @warnIfNoDataLoaded
    @updateDisplay
    def sort(self):
        modelops.sortTrajs(self.data)

    @warnIfNoDataLoaded
    @updateDisplay
    def slowSort(self):
        modelops.sortTrajsSlow(self.data)

    @warnIfNoDataLoaded
    def plotContinuity(self):
        analysis.plotContinuity(self.data)

    @warnIfNoDataLoaded
    def plotEnergy(self):
        analysis.plotEnergy(self.data)

    @warnIfNoDataLoaded
    def plotSpeed(self):
        analysis.plotSpeed(self.data)

    @warnIfNoDataLoaded
    def plotLengthHistogram(self):
        analysis.plotLengthHistAll(self.data)

    @warnIfNoDataLoaded
    @updateStatus
    def fillSmallGaps(self):
        progress = self.mkProgress("Filling small gaps...")
        modelops.fillGaps(self.data, 0.2, 0.2, progress)

    @warnIfNoDataLoaded
    @updateStatus
    def fillAllGaps(self):
        progress = self.mkProgress("Filling all gaps...")
        modelops.fillGaps(self.data, self.data.numFrames / self.data.framerate, 0.2, progress)

    @warnIfNoDataLoaded
    @updateNumFrames
    @updateStatus
    def cutRight(self):
        self.model.cutRight()

    @warnIfNoDataLoaded
    @updateNumFrames
    @updateStatus
    def cutLeft(self):
        self.model.cutLeft()

    @warnIfNoDataLoaded
    @updateStatus
    def removeUnassigned(self):
        self.model.deleteUnassigned()

    # End actions #

    def exitNoSave(self):
        if self.data is not None:
            changed = self.data.changed
            if changed:
                reply = QtGui.QMessageBox.question(
                    self, 'Message', "Are you sure you want to lose unsaved data?", 
                    QtGui.QMessageBox.Yes | QtGui.QMessageBox.Cancel,
                    QtGui.QMessageBox.Cancel)
                if reply == QtGui.QMessageBox.Cancel:
                    return
            self.data.changed = False
        self.close()
        if self.data is not None:
            self.data.changed = changed

    def closeEvent(self, event):
        self.saveFile()
        event.accept()

    def keyPressEvent(self, e):
        self.gl.keyPressEvent(e)

    @updateDisplay
    @updateStatus
    def setFrame(self, numFrame):
        self.model.frame = numFrame

    def updateFrame(self):
        self.transport.setFrame(self.model.frame)

    @ifDataLoaded
    def updateStatus(self):
        msg = "Frame: %d | traj: %d of %d | " % \
            (self.model.frame, self.model.traj, self.data.numTrajs)
        msg += self.model.humanReadPart()
        traj = self.model.currentTraj()
        frame = self.model.frame
        if (traj.hasFrame(frame)):
            point = traj.getFrame(frame)
            msg += " at (%.0f, %.0f, %.0f)" % tuple(point)
            self.statusBar().showMessage(msg)
        else:
            msg += " (not in this frame)"
            self.statusBar().showMessage(msg)

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    cd = CuteDotsMainWindow()

    def excepthook(extype, exval, tb):
        message = "An error has occured. Technical details:\n\n" + str(exval)
        traceback.print_tb(tb)
        QtGui.QMessageBox.critical(cd, "Error", message + str(exval))
    sys.excepthook = excepthook

    cd.show()
    sys.exit(app.exec_())    
