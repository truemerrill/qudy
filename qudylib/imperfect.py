# IMPERFECT.PY
#
# A class for imperfect propagators
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

from quantop import *
from propagator import *
import error, control, copy


__all__ = ['imperfect']


class imperfect( propagator ):
    """
    A class for imperfect quantum dynamical propagators.  Similar to
    the propagator class, however now we associate an error model that
    distorts the input control functions.
    
    """
    
    def __init__(self, *args, **keyword_args):
        
        # Pull out last argument.  It should be an instance of the
        # error class.  If it instead a string, interpret the string
        # as a name of an error function.
        
        err = args[ len(args) - 1 ]
        if isinstance( err, error.error ):
            self.error = err
            
        elif isinstance( err, str ):
            self.error = error.error( err )
            
        else:
            raise SyntaxError("Final argument must be an error object.")
        
        # For the remaining arguments, use propagator's __init__
        # function to produce the relevant methods and objects.
        
        propagator.__init__( self, *args[0:len(args)-1], **keyword_args )

        # Save an ideal set of controls
        self.ideal_control = copy.copy( self.control )
        
        # Update the error
        self.update_error()
        
        
    def __repr__(self):
        """
        Function to display imperfect objects when called on the
        command line.
        """
        string = str( self.dimension ) + '-D imperfect propagator on t = ( ' + \
            str( self.timemin ) + ' , ' + str( self.timemax ) + ' )' 
        return string
    
    
    def __mul__(self, target):
        """
        Function to multiply two imperfect objects.  For propagators
        for the same bilinear control system, i.e. their
        `self.hamiltonians` entries match, multiplication involves
        appending the control functions to make a larger control.
        """
        
        # Use propagator __mul__ method, but make the output an
        # imperfect instance and make the error match self.
        
        U = propagator.__mul__(self, target)
        V = imperfect( U.ideal_control , self.error )
        return V
    
    
    def update_error(self, *args):
        """
        Function to update the error model.  Using the new error model, the
        control object is replaced with it's distorted counterpart, ie.
        
        # Insert math stuff here
        
        **Forms:**
        
           * `update_error()`
           * `update_error( err )`
           * `update_error( model )`
           * `update_error( model, error_parameters )`
           
        **Args:**
        
           * *err* : An error object.
           * *model* : A string representing the name of the error model
           * *error_parameters* : An ordered list of parameters required
             by the error model.  When the parameters are not provided,
             the model uses a default set. Check the error model's
             specific documentation for more details
             
        When no arguments are supplied, the method uses the error instance
        currently mapped to imperfect.
        """
        
        if not len(args) == 0:
        
            if isinstance( args[0] , error.error ):
                # Input is a new error model, replace
                self.error = args[0]
                
            else:
                # Input is parameters for a new error model
                err = error.error( *args )
                self.error = error
                
            # Recalculate with new error model
            self.update_error()
            
        else:
            
            # Replace self.control with distorted controls
            #ideal = control.control( self.ideal_control, self.times )
            distorted = self.error( self.ideal_control )
            self.control = distorted
