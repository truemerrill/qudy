# AMPLITUDE.PY
#
# Amplitude error model
# 
# Copyright (C) 2011, 2012 True Merrill
# Georgia Institute of Technology
# School of Chemistry and Biochemistry
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from ..control import control

def call( ctrl, error_parameters, **keyword_args ):
    """
    Method for amplitude error model.
    """
    
    # Enforce defaults
    if error_parameters == None:
        error_parameters = default_parameters()
      
    # Perhaps the user was lasy any input the error parameters as a
    # single element rather than a subscriptable list.  Fix.
    try:
        epsilon = error_parameters[0]
    except TypeError:
        epsilon = error_parameters
    
    # Update amplitude of control values.
    arr = (1.0 + epsilon) * ctrl.control
    t = ctrl.times

    # Return modified control.
    return control(arr,t)


def default_parameters():
    """
    Default parameters.
    """
    
    return [0.1]


def repr( error_parameters ):
    """
    Function to display amplitude objects when called on the command line
    """
    
    string = "amplitude error: \n" + \
             "    epsilon:\t%.2E" %( error_parameters[0] )
    return string
