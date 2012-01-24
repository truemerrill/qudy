# PROPAGATOR.PY
#
# A class for solutions to bilinear control problems
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
from integration import *
from routines import *


class propagator:
    """
    A class for quantum dynamical propagators.  The evolution is
    specified by a bilinear control system, and optionally, a set of
    Linblad decay channels.
    
    Inputs may be an instance of the control class and also lists of N
    x N matrices, where N is the dimensionality of the Hilbert space.
    If no Lindblad channels are specified, then the evolution is
    assumed to be unitary.  If no Hamiltonians are specified and the
    dimensionality of the control system is n(n+1), then we assume the
    dynamical Lie algebra is SU(2^n) and use a product-operator basis
    for the
    
    **Forms:**
    
       * `propagator(ctrl, hamiltonians)`
       * `propagator(ctrl, hamiltonians, solution = 'method')`
       * `propagator(ctrl, hamiltonians, solution = 'method', order = n)`
    
    **Args:**
      
       * *ctrl* : An instance of the control class.  Contains time
          information as well as k-many control functions.
       * *hamiltonians* :  A list or array of k-many Hamiltonians.  
          The Hamiltonians must be square matrices of the same 
          dimensionality.  It is recommended to use the operator 
          class.
    
    **Optional keywords:**
    
       * solution = 'method' : Method used to numerically solve the 
         conrol problem.  'method' may be one of the following strings.
         
            1. 'trotter' : Trotter method, using each time interval as
               a time slice.
            2. 'dyson' : Dyson series.
            3. 'magnus' : Magnus expansion.
            4. 'lindblad' : Lindblad master equation.
            
       * order = n : Integer valued order of pertubaton theory.  Used in 
         the Dyson and Magnus methods
       * integration_method = 'method' : Integration technique used in 
         the Dyson and Magnus methods. 'method' may be one of the following
         stings.
                            
            1. 'trapz' :    trapezoidal rule
            2. 'cumtrapz' : cumulative trapezoidal
            3. 'romb' :     Romberg integration
            4. 'simps' :    Simpson's rule
            
    **Returns:**
    
       * :math:`U`, matrix solution to the control equation 
         :math:`\dot{U}(t) = \sum_\mu u_\mu(t) H_\mu U(t)`.
    """
    def __init__(self, *args, **keyword_args):
        
        # If there is one argument, then the user only gave a control
        # function.  Assume a qubit system, i.e. SU(2^n).
        if len( args ) ==  1:
            ctrl = args[0]
            
            try:
                dim = ctrl.dimension
                # Invert Lie dimensionality to get number qubits
                n = log( dim + 1 ) / log( 4 )
                H = product_operator( n )
                
            except AttributeError:
                raise TypeError('ctrl must be a member of the control class.')
            
            self.Lie_algebra = 'su(' + str(2**n) + ')'
            self.number_qubits = n
            
        # If there are two arguments, then the user gave both a
        # control function and also a set of Hamiltonians.
        elif len( args ) == 2:
            ctrl = args[0]
            H = args[1]
            
            try:
                dim = ctrl.dimension
                dimH = len(H)
                
            except AttributeError:
                raise TypeError('ctrl must be a member of the control class.')
            
            # Check that ctrl and H have the same dimension
            if not len(ctrl) == len(H):
                raise ValueError('Dimension of ctrl %i and dimension of ' + \
                      'hamiltonians %i do not match.' %(int(dim),int(dimH)))
        
        # If there are three arguments, then the user gave a control
        # function, Hamiltonians and also a set of Lindblad operators.
        elif len( args ) == 3:
            ctrl = args[0]
            H = args[1]
            L = args[2]
            
            try:
                dim = ctrl.dimension
                dimH = len(H)
                
            except AttributeError:
                raise TypeError('ctrl must be a member of the control class.')
            
            # Check that ctrl and H have the same dimension
            if not len(ctrl) == len(H):
                raise ValueError('Dimension of ctrl %i and dimension of ' + \
                      'hamiltonians %i do not match.' %(int(dim),int(dimH)))
            
            # Add in Lindblad operators and force solution method to
            # be Lindblad.
            self.lindblad = L
            
           
        # Did not understand inputs, too many!
        else:
            raise SyntaxError('Too many inputs.')
        
        self.control = ctrl
        self.hamiltonians = H
        self.dimension = self.control.dimension
        self.number_controls = self.control.number_controls
        self.timemin = self.control.timemin
        self.timemax = self.control.timemax
        
        # Parse through keyword arguments.  Sets default solution
        # method.  Other keywords that are not understood will be
        # quietly ignored.
        if keyword_args.has_key( 'solution' ):
            
            method = keyword_args['solution']
            valid_inputs = ['trotter', 'dyson', 'magnus', 'lindblad']
            
            if method in valid_inputs:
                self.solution_method = method
                
            else:
                # Keyword was not valid.  Check to see if already
                # defined.
                
                # Return to default and warn user
                self.solution_method = 't=[rotter'
                warn('Solution method not understood, defaulting to \'trotter\'.')
        else:
            # set default solution method
            self.solution_method = 'trotter'
            
        # Set order of pertubation theory
        default_order = 4
        if keyword_args.has_key( 'order' ):
            
            order = keword_args['order']
            
            try:
                if (order % 1.0) == 0:
                    self.order = order
            
                else:
                    # Order was not an integer, revert to defaults
                    warn('Pertubation theory order must be integer valued.')
                    self.order = default_order
                    
            except TypeError:
                # Order was not numeric, revert to defaults
                warn('Pertubation theory order must be integer valued.')
                self.order = default_order
                
        else:
            # Set default order
            self.order = default_order
            
                       
    
    def __repr__(self):
        """
        Function to display propagator objects when called on the
        command line.
        """
        string = str( self.dimension ) + '-D propagator on t = ( ' + \
                 str( self.timemin ) + ' , ' + str( self.timemax ) + ' )'
        return propagator( )


    def __mul__(self, target):
        """
        Function to multiply two propagator objects.  For propagators
        for the same bilinear control system, i.e. their
        `self.hamiltonians` entries match, multiplication involves
        appending the control functions to make a larger control.
        """
        
        def H_check( h1, h2 ):
            # todo : fix so that hamiltonians are checked.
            return h1 == h2
        
        try:
            if H_check( self.hamiltonains , target.hamiltonians ):
                
                # Hamiltonains match, append controls together
                ctrl1 = target.control
                ctrl2 = self.control
                ctrl = vstack( ctrl1, ctrl2 )
                
                # Adjust time vectors.  The second control should
                # start just as the first one finishes.
                t1 = target.times
                t2 = self.times
                
                # To avoid the control function from having two values
                # at the same instant in time, we add an infinitesimal
                # delay.
                dt_min = 1E-10 * target.timemax
                
                t2 = (t2 - self.timemin) + target.timemax + dt_min
                t = vstack( t1, t2 )
                
                # Stack control and time vectors.
                ARR = hstack( ctrl, t )
                
            else:
                raise ValueError('Hamiltonians do not match. ' +\
                      'Multiplication is ill-defined.')
        
        except AttributeError:
            raise TypeError('Multiplication is only defined between ' +\
                      'propagator objects.')
        
        return propagator( control(ARR) , self.hamiltonains )
    
    
    def solve(self, method = 'trotter'):
        """
        
        """
        
        if method == 'trotter':
            U = trotter( self.control, self.hamiltonians )
            
        elif method == 'dyson':
            U = dyson( self.control, self.hamiltonains, \
                       self.order )
            
        elif method == 'magnus':
            U = magnus( self.control, self.hamiltonains, \
                        self.order )
            
        elif method == 'lindblad':
            U = lindblad( self.control, self.hamiltonians, \
                          self.lindblad )
            
        else:
            raise ValueError('Method %s was not understood.' %method)
        
        return U



