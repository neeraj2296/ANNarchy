"""
    
    ANNarchyEditor.py
    
    This file is part of ANNarchy.
    
    Copyright (C) 2013-2016  Julien Vitay <julien.vitay@gmail.com>,
    Helge Uelo Dinkelbach <helge.dinkelbach@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ANNarchy is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
"""
from PyQt4.QtGui import QApplication, QMainWindow
from PyQt4.QtCore import pyqtSignal, pyqtSlot, QObject, SIGNAL

import sys
import ANNarchy4
from Repository import Repository
from VisWidget import VisualizerWidget

# generated by pyuic4
from ui_ANNarchyEditor import Ui_ANNarchyEditor 
            
class EditorMainWindow(QMainWindow):
    """
    Main class of editor subpackage. It represents the main widget and instantiate
    an UI class Ui_ANNarchyEditor created by the pyuic4 tool.
    """
    initialize = pyqtSignal()
    signal_change_grid_from_general = pyqtSignal(int, int)
    signal_net_editor_to_pop_view = pyqtSignal(int, int)
    
    def __init__(self, func):
        """
        Constructor initializes first all widgets in the hierarchy and then 
        emits the ``initialize`` signal to call the corresponding function in 
        all subwidgets.
        
        Parameter:
            
        * *func*: user defined function, defining the simulation loop.
        """
        super(QMainWindow, self).__init__()

        ANNarchy4.Global._visualizer = self
        
        self._rep = Repository()
        self._rep.load()
        
        self.setWindowTitle('ANNarchy4.1 ultimate editor') 
        
        self._ui = Ui_ANNarchyEditor()
        self._ui.setupUi(self)
        self._func = func

        self._connect_signals()

        #
        # initialize all child widgets
        self.emit(SIGNAL("set_repository(PyQt_PyObject)"), self._rep)
        self.initialize.emit()
        
        self._visualizer = self.findChild(VisualizerWidget, 'visualizer')

    def _connect_signals(self):
        """
        Signals for standard types and Qt objects are working fine by the editor tool. For own classes
        etc. we need the old style of SIGNAL emitation.
        
        .. hint: please note, that the receiving functions are not decorated by @pyqtSlot
        """
        QObject.connect(self, SIGNAL("set_repository(PyQt_PyObject)"), self._ui.objects.set_repository)
        QObject.connect(self, SIGNAL("set_repository(PyQt_PyObject)"), self._ui.neur_general.set_repository)
        QObject.connect(self, SIGNAL("set_repository(PyQt_PyObject)"), self._ui.syn_general.set_repository)
        QObject.connect(self, SIGNAL("set_repository(PyQt_PyObject)"), self._ui.net_select.set_repository)
        QObject.connect(self, SIGNAL("set_repository(PyQt_PyObject)"), self._ui.editor.set_repository)        
        QObject.connect(self, SIGNAL("set_repository(PyQt_PyObject)"), self._ui.pop_view.set_repository)
        QObject.connect(self, SIGNAL("set_repository(PyQt_PyObject)"), self._ui.complete.set_repository)
        
    def set_data(self, x, y, x_data, y_data):
        """
        shortcut function for tests, to be deleted later
        """
        self._visualizer.set_data(x, y, x_data, y_data)

    def render_data(self):
        self._visualizer.render_data()

    @pyqtSlot()
    def compile_and_run(self):
        """
        pyqtSlot function.
        
        Compiles the ANNarchy library and runs the simulation loop.
        
        Parameter:
        
        * *none*
        
        Emitted by:
        
        * Mainwindow compile button pressed. 
        """
        ANNarchy4.compile()
        
        self._func()

    def closeEvent(self, event):
        """
        Overloaded PyQt4.QtGui.QMainWindow function.
        
        The editor was closed, so we store the repository data and quit.
        """
        self._rep.save()
        
        return QMainWindow.closeEvent(self, event)
        
    def resizeEvent(self, resize_event):
        """
        Overloaded PyQt4.QtGui.QMainWindow function.

        The window was resized and new border sizes for the TabWidget are determined.
        The TabWidget should take at least 90 percent of width and 80 percent of window
        height.
        
        Parameter:
            
        * *resize_event*: not needed by this function, but piped to parent class
        """
        new_border_w = int (self.width() * 0.80)
        new_border_h = int (self.height() * 0.90)
        
        self._ui.splitter.setMinimumSize( new_border_w, new_border_h )
        self._ui.splitter.setMaximumSize( new_border_w+1, new_border_w )

        return QMainWindow.resizeEvent(self, resize_event)
            
class ANNarchyEditor(object):
    """
    Entry point for the editor module. Instantiates the EditorMainwindow
    and after this enter the QT main loop.
    """
    def __init__(self, func):
        app = QApplication(sys.argv)
         
        win = EditorMainWindow(func)
        
        win.show()
         
        sys.exit(app.exec_())
        