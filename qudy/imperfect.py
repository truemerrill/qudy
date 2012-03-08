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
import error, control


__all__ = ['imperfect','imperfect_rotation','M']


class imperfect( propagator ):
    """
    A class for imperfect quantum dynamical propagators.  Similar to
    the propagator class, however now we associate an error model that
    distorts the input control functions.
    
    **Forms:**
    
       * ``imperfect(ctrl, err)``
       * ``imperfect(ctrl, 'error model')``
       * ``imperfect(ctrl, hamiltonians, err)``

    **Args:**
      
       * *ctrl* : An instance of the control class.  Contains time
          information as well as k-many control functions.
       * *'error model'* a string corresponding to the name of an
         error model.
       * *hamiltonians* :  A list or array of k-many Hamiltonians.  
         The Hamiltonians must be square matrices of the same 
         dimensionality.  It is recommended to use the operator 
         class.
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
        self.ideal_control = self.control.copy()
        
        # Update the error
        self.update_error()
        
        
    def __repr__(self):
        """
        Function to display imperfect objects when called on the
        command line.
        """
        string = str( self.dimension ) + '-D imperfect propagator on t = ( ' + \
            str( self.timemin() ) + ' , ' + str( self.timemax() ) + ' )' 
        return string
    
    
    def __mul__(self, target):
        """
        Function to multiply two imperfect objects.  For propagators
        for the same bilinear control system, i.e. their
        ``self.hamiltonians`` entries match, multiplication involves
        appending the control functions to make a larger control.
        """
        
        # Use propagator __mul__ method, but make the output an
        # imperfect instance and make the error match self.
        U1 = propagator( self.ideal_control )
        U2 = propagator( target.ideal_control )
        U = propagator.__mul__(U1,U2)
        
        V = imperfect( U.ideal_control , self.hamiltonians, self.error )
        return V
    
    
    def copy(self):
        """
        Create an independent copy of self in memory.
        """
        ctrl = self.ideal_control.copy()
        hamiltonians = self.hamiltonians[:]
        error = self.error.copy()

        c = imperfect( ctrl, hamiltonians, error )
        c.solution_method = str( copy( self.solution_method ) )
        c.order = int( copy( self.order ) )
        
        return c
    
    
    def inverse(self):
        """
        Calculates the propagator inverse of self.  Note that unlike
        in matrix methods, which would require diagonalizing the
        propagator, here we compute the inverse by a reording of the
        control functions (much faster).
        """
        
        # Calculate inverse control
        ctrl = self.ideal_control.control.copy()
        inv_ctrl = - flipud( ctrl )
        
        # Create copy of self
        c = self.copy()
        
        # Replace c.control with inverse controls
        c.ideal_control.control = inv_ctrl
        c.update_error()

        return c

        
    def update_error(self, *args):
        """
        Function to update the error model.  Using the new error model, the
        control object is replaced with it's distorted counterpart, ie.
        
        # Insert math stuff here
        
        **Forms:**
        
           * ``update_error()``
           * ``update_error( err )``
           * ``update_error( model )``
           * ``update_error( model, error_parameters )``
           
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


    def interaction_frame( self ):
        """
        Transform into interaction frame moving with the
        ideal_control.  The transformation into the interaction
        picture is an operation on the basis Hamiltonians, which then
        can be mapped onto an image of the control functions
        themselves.
        """
        return NotImplemented


def imperfect_rotation( *args, **keyword_args ):
    """
    A function to form propagators that represent rotations in SU(2).
    These may be multiplied to produce a pulse sequence.
    
    **Forms:**
    
       * ``imperfect_rotation( axis, err )``
       * ``imperfect_rotation( theta, phi, err )``
       * ``imperfect_rotation( theta, axis, err )``
       
    **Args:**
      
       * *axis* : A three-element list, tuple or array representing
         components of a Bloch vector.  If axis is the sole input,
         then the rotation angle is interpreted to be the length of
         the axis vector.
       * *theta* : A rotation angle.
       * *phi* : A field phase.  The interaction frame Hamiltonian for
         this phase is proportional to:math:`H = \\cos \\phi X + 
         \\sin \\phi Y`.

    **Optional keywords:**
          
       * See the propagator class for keywords.
     
    **Returns:**
    
       * :math:`V = M(\\theta,\\phi)`,
    """

    # Pull out last argument.  It should be an instance of the
    # error class.  If it instead a string, interpret the string
    # as a name of an error function.
        
    err = args[ len(args) - 1 ]

    # Construct ideal rotation, steal control to make an imperfect
    # propagator.
    U = rotation( *args[0:(len(args)-1)], **keyword_args )
    ctrl = U.control

    return imperfect( ctrl, err )


def M( *args, **keyword_args ):
    """
    A convenient shortcut for imperfect_rotation().  See
    imperfect_rotation().
    """
    return imperfect_rotation( *args, **keyword_args )
