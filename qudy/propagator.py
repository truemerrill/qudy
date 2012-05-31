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
import control
import integration
import routines
import imperfect


__all__ = ['propagator','rotation','R']

[Hx,Hy,Hz] = routines.product_operator(1)

class propagator:
    """
    A class for quantum dynamical propagators.  The evolution is
    specified by a bilinear control system, and optionally, a set of
    Linblad decay channels.
    
    Inputs may be an instance of the control class and also lists of N
    x N matrices, where N is the dimensionality of the Hilbert space.
    If no Lindblad channels are specified, then the evolution is
    assumed to be unitary.  If no Hamiltonians are specified and the
    dimensionality of the control system is n(n-1), then we assume the
    dynamical Lie algebra is SU(2^n) and use a product-operator basis
    for the
    
    **Forms:**
    
       * ``propagator(ctrl)``
       * ``propagator(ctrl, hamiltonians)``
       * ``propagator(ctrl, hamiltonians, solution = 'method')``
       * ``propagator(ctrl, hamiltonians, solution = 'method', order = n)``
    
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
    
       * :math:`U`, propagator representing the solution to the
         control equation :math:`\dot{U}(t) = \sum_\mu u_\mu(t) H_\mu
         U(t)`.
    """
    def __init__(self, *args, **keyword_args):
        
        # If there is one argument, then the user only gave a control
        # function.  Assume a qubit system, i.e. SU(2^n).
        if len(args) == 1:
            self.ideal_control = args[0]
            dim = self.ideal_control.number_controls
                
            # Convert Lie dimensonality into number of qubits.
            n = log( dim + 1 ) / log( 4 )
            self.number_qubits = n
            self.lie_algebra = "su(%i)" %(int(2**n))
            if n == 1:
                self.hamiltonians = [Hx,Hy,Hz]
            else:
                self.hamiltonians = routines.product_operator( n )
            
        # If there are two arguments, then the user gave both a set of
        # controls and also a set of Hamiltonians.
        elif len(args) == 2:
            self.ideal_control = args[0]
            self.hamiltonians = args[1]
            
        # If there are three arguments, then the user gave controls,
        # Hamiltonians and Lindblad operators.
        elif len(args) == 3:
            self.ideal_control = args[0]
            self.hamiltonians = args[1]
            self.lindblad = args[2]
            
            # Since Lindblad operators were specified, we should use
            # Master equation methods.  Calculate superoperator
            # representations of relevant operators.
            
        if not len( self.ideal_control ) == len( self.hamiltonians ):
            raise ValueError("Control dimension mismatch.  Controls and " + \
                             "Hamiltonians must be of the same length.")
            
        # Carry over several constants from controls
        self.dimension  =  self.ideal_control.dimension
        self.times      =  self.ideal_control.times
        self.timemin    =  self.ideal_control.timemin
        self.timemax    =  self.ideal_control.timemax
        
        # Create an control.  For propagator instances, these are
        # identical to ideal_control, however for imperfect instances
        # these differ.
        self.control = self.ideal_control.copy()
        
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
                self.solution_method = 'trotter'
                warn('Solution method not understood, defaulting to \'trotter\'.')
        else:
            # set default solution method
            self.solution_method = 'trotter'
            
        # Set order of pertubation theory
        default_order = 4
        if keyword_args.has_key( 'order' ):
            
            order = keyword_args['order']
            
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
                 str( self.timemin() ) + ' , ' + str( self.timemax() ) + ' )'
        return string


    def __mul__(self, target):
        """
        Function to multiply two propagator objects.  For propagators
        for the same bilinear control system, i.e. their
        `self.hamiltonians` entries match, multiplication involves
        appending the control functions to make a larger control.
        """
        
        def H_check( h1, h2 ):
            # Cleverly uses a method in the operator class to
            # determine if Hamiltonians are identical.
            return h1 == h2
        
        try:
            if H_check( self.hamiltonians , target.hamiltonians ):
                
                # Hamiltonians match, append controls together
                ctrl1 = target.ideal_control.control
                ctrl2 = self.ideal_control.control
                ctrl0 = vstack( (ctrl1, ctrl2) )
                
                # Adjust time vectors.  The second control should
                # start just as the first one finishes.
                t1 = target.times
                t2 = self.times
                
                # To avoid the control function from having two values
                # at the same instant in time, we add an infinitesimal
                # delay.
                dt_min = 1E-10 * target.timemax()
                
                t2 = (t2 - self.timemin() ) + target.timemax() + dt_min
                t = vstack( (t1, t2) )
                
                # Stack control and time vectors.
                ARR = hstack( (ctrl0, t) )
                
            else:
                raise ValueError('Hamiltonians do not match. ' +\
                      'Multiplication is ill-defined.')
        
        except AttributeError:
            raise TypeError('Multiplication is only defined between ' +\
                      'propagator objects.')

        if isinstance(target, imperfect.imperfect):
            # Check if target is an imperfect propagator.  If it is,
            # we should default to returning an imperfect propagator.
            return imperfect.imperfect( control.control(ARR) , \
                               self.hamiltonians, target.error )
        
        else:
            # Must be a normal propagator.
            return propagator( control.control(ARR) , self.hamiltonians )

    
    def __call__(self, time = None):
        """
        
        """
        
        if time == None:
            return self.solve()
        
        elif time > self.timemax() or time < self.timemin():
            raise ValueError('Time is not within interval bounds' + \
                  ' ( %.2E , %.2E ).' %(self.timemin(), self.timemax() ))
        
        try:
            t = min( self.control.times[ self.control.times - time > 0 ] )
        except ValueError:
            # We must be at the higher time limit
            t = time
            
        index = self.times.tolist().index( t )
        
        # Cut controls over interval (timemin, time)
        arr = self.control.control[ 0:index+1 , : ]
        times = self.times[ 0:index+1 ]
        ARR = hstack( (arr,times) )
        
        # Form new propagator
        U = propagator( control.control( ARR ) )
        
        # Solve propagator
        return U.solve()
    
    
    def copy(self):
        """
        Creates an independent copy of self in memory.
        """
        
        ctrl = self.control.copy()
        hamiltonians = self.hamiltonians[:]
        c = propagator( ctrl, hamiltonians )
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
        ctrl = self.control.inverse()
        
        # Create copy of self
        c = self.copy()
        
        # Replace c.control with inverse controls
        c.control = ctrl
        c.ideal_control = ctrl
        
        return c
        
    
    def solve(self, method = 'trotter'):
        """
        
        """
        
        if method == 'trotter':
            U = integration.trotter( self.control, \
                self.hamiltonians )
            
        elif method == 'dyson':
            U = integration.dyson( self.control, \
                self.hamiltonians, self.order )
            
        elif method == 'magnus':
            U = integration.magnus( self.control, \
                self.hamiltonians, self.order )
            
        elif method == 'lindblad':
            U = integration.lindblad( self.control, \
                self.hamiltonians, self.lindblad )
            
        else:
            raise ValueError('Method %s was not understood.' %method)
        
        return U


    def components(self, *args):
        """
        Calculates components of generator on the Lie algebra.
        """

        # Skew-symmetrize current Hamiltonians
        H = []
        for index in range(len( self.hamiltonians )):
            H.append( -1j * self.hamiltonians[index] )

        # Generate an orthanormal set of skew-symmetrized Hamiltonians
        # that span the Lie algebra.
        A = routines.generate_algebra( self.hamiltonians )

        # Decompose log(U) into A basis
        U = self.solve()
        logU = operator( logm(U) )
        
        a = zeros( len(A) )
        for v in range( len(A) ):
            a[v] = routines.inner_product( logU, A[v] )

        # Calculate covectors in original oblique basis
        print a
        u = zeros( len(H) )
        for mu in range(len(H)):
            
            for v in range( len(A) ):
                u[v] = u[v] + a[v] / routines.inner_product(A[v],H[mu])

        return u[v]
    

def rotation( *args, **keyword_args ):
    """
    A function to form propagators that represent rotations in SU(2).
    These may be multiplied to produce a pulse sequence.
    
     **Forms:**
    
       * ``rotation( axis )``
       * ``rotation( theta, phi )``
       * ``rotation( theta, axis )``
       
    **Args:**
    
       * *axis* : A three-element list, tuple or array representing
         components of a Bloch vector.  If axis is the sole input,
         then the rotation angle is interpreted to be the length of
         the axis vector.
       * *theta* : A rotation angle.
       * *phi* : A field phase.  The interaction frame Hamiltonian for
         this phase is proportional to :math:`H = \\cos \\phi X + 
         \\sin \\phi Y`.
       
    **Optional keywords:**
    
    See the propagator class for keywords.
     
    **Returns:**
    
       * :math:`U = R(\\theta,\\phi)`,
    """
    
    # Check input arguments
    if len( args ) == 1:
        
        # Single argument, so must have been the rotation axis. Check
        # if argument is a 3-vector.
        
        arg = args[0]
        try:
            
            if arg.__len__() == 3:
                 theta = sqrt( arg[0]**2 + arg[1]**2 + arg[2]**2 )
                 vec  = array([ arg[0], arg[1], arg[2] ]) / theta
                 
            else:
                raise SyntaxError("Axis must be a 3-vector.")
            
        except AttributeError:
            raise SyntaxError("Axis must be a 3-vector.")
        
    elif len( args ) == 2:
        
        # Either an angle and an axis or an angle and a phase.
        # Distinguish between the two and construct a 3-vector.
        
        theta = args[0]
        arg = args[1]
        try:
            
            if arg.__len__() == 3:
                # Normalize the axis
                norm = sqrt( arg[0]**2 + arg[1]**2 + arg[2]**2 )

                if norm == 0:
                    raise ValueError("Axis must not be a null vector.")

                vec  = array([ arg[0], arg[1], arg[2] ]) / (norm)
                
            else:
                raise SyntaxError("Axis must be a 3-vector.")
            
        except AttributeError:
            # Must have specified a phase instead of an axis.  Produce
            # one for the user.
            
            phi = arg
            vec = array([cos(phi), sin(phi), 0])
            
    else:
        raise SyntaxError("Improper number of arguments.")
    
    # Construct the control function and the propagator.  We assume
    # square pulse propagators (i.e. constant controls over the
    # interval).
    
    # Check for negative rotation angles.
    if theta < 0:
        vec = - vec
        theta = - theta
    
    ctrl = array([[vec[0], vec[1], vec[2], 0    ],
                  [vec[0], vec[1], vec[2], theta]])
    
    u = control.control(ctrl)
    return  propagator( u, [Hx,Hy,Hz], **keyword_args )


def R( *args, **keyword_args ):
    """
    A convenient shortcut for rotation().  See rotation().
    """
    return rotation( *args, **keyword_args )
    
    
                
        
            
