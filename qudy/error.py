# ERROR.PY
#
# A class for quantum control functions
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
from control import control

__all__ = ['error']

class error:
    """
    class for quantum control error models.
    
    An error model describes deformations to the control functions
    caused by noise or experimental imperfections.  Each error model
    is implemented as an instance of the error class, and is
    initialized by specifying the model name and a set of error
    parameters.  
    
    **Forms:**
    
       * ``error( model )``
       * ``error( model, error_parameters )``
       
    **Args:**
    
       * *model* : A string representing the name of the error model
       * *error_parameters* : An ordered list of parameters required
         by the error model.  When the parameters are not provided,
         the model uses a default set. Check the error model's
         specific documentation for more details.
    
    Once constructed, an error model acts as a functional: taking an
    input control and returning an output (distorted) control.  
    Consider the following code.
    
    .. code-block:: python
       
       # Construct a control instance
       u1 = lambda x: sin(x)
       u2 = lambda x: cos(x)
       t = arange( 0, 2*pi, pi/50.0 )
       u = control( u1, u2, t )
       
       # Contruct an error instance
       err = error( 'amplitude', [0.1] )
       
       # Return distorted controls
       v = err( u )
    
    In certain cases, such as when the model uses a master equation to
    compute the evolution, the error instance will return both a set
    of distorted controls and a set of Lindblad operators.  These may
    be used as inputs to a propagator instance.  Consider the
    following code.
    
    .. code-block:: python
       
       # Construct an error instance
       err = error( 'depolarizing' )
       
       # Return distorted controls
       [v, Lindblad] = err(u)
       
    """
    def __init__( self, *args, **keyword_args ):
        """
        Initialize an error instance.
        """
        
        # First argument should always be the model name
        model_name = args[0]
        
        try:
            exec('import error_models.%s as model' %(model_name))
            
        except ImportError:
            raise ImportError('Could not find error model %s.' %(model_name))
        
        # Save the model and model_name
        self.model = model
        self.model_name = model_name
        
        # Second argument represents the error parameters.  If the
        # second argument does not exist, then pull the default
        # parameters from model.
        if len( args ) == 1:
            error_parameters = self.model.default_parameters()
            
        else:
            error_parameters = args[1]
            
            if not isinstance(error_parameters, list):

                try:
                    error_parameters = list(error_parameters)
                except TypeError:
                    error_parameters = [ float(error_parameters) ]
            
        # Save the error parameters
        self.default_parameters = self.model.default_parameters()
        self.error_parameters = error_parameters
        
        
    def __call__( self, ctrl ):
        """
        Implements functional behavior.  Inputs a control instance and
        returns a distorted control function.
        """
        return self.model.call( ctrl, self.error_parameters )
    
    
    def __repr__( self ):
        """
        Function to display error objects when called on the command line.
        """
        try:
            return self.model.repr( self.error_parameters )
        
        except AttributeError:
            
            # Error model must have left off a repr routine.  Use a
            # generic method.
            
            num_epsilons = len( self.error_parameters )
            string = "unspecified error: "
            
            for index in range(num_epsilons):
                s = "\n    epsilon%i:\t%.2E" %( index+1 , self.error_parameters[index] )
                string = string + s
                
            return string
        
        
    def copy(self):
        """
        Function to make copy of self in memory.
        """
        return error( self.model_name, self.error_parameters )
