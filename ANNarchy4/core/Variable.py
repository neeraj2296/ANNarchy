"""

    Variable.py
    
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
class Variable(object):
    """
    Variable representation in ANNarchy.
    """
    def __init__(self,  **keyValueArgs):
        """
        Set of key-value pairs defining this variable.
        
        Possible keys are:
        
            * *eq*: variable equation
            * *init*: initialization
            * *min*: minimal value for this variable
            * *max*: maximal value for this variable            
            
        other keys will be rejected.
        """
        self.eq = None
        self.init = None
        self.min = None
        self.max = None
        
        for key in keyValueArgs:

            if key == 'eq':
                self.eq = keyValueArgs[key]
            elif key=='init':
                self.init = keyValueArgs[key]
            elif key=='min':
                self.min = keyValueArgs[key]
            elif key=='max':
                self.max = keyValueArgs[key]
            else:
                print 'unknown key: '+key
         