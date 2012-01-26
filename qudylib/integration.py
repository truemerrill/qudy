# INTEGRATION.PY
# 
# Integration subroutines for quantum control problems
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
from scipy.integrate import trapz, cumtrapz, simps, romb


def integrate( ctrl, hamiltonians, method = 'trapz' ):
    """
    Integrates a bilinear control system consisting of a control and a
    set of Hamiltonians.  Several numerical integration routines are
    avaliable.  Wrapper for the scipy.integrate backend.
    
    **Forms:**
    
        * ``integrate( ctrl, hamiltonians )``
        * ``integrate( ctrl, hamiltonians, method = 'method' )``
        
    **Args:**
    
        * *ctrl* : An instance of the control class.  Contains time
          information as well as k-many control functions.
        * *hamiltonians* :  A list or array of k-many Hamiltonians.  
          The Hamiltonians must be square matrices of the same 
          dimensionality.  It is recommended to use the operator 
          class.
                 
    **Optional keys:**
    
        * method = 'method': Chooses an integration method.  The 
          method may be one of the following options.
                 
             1.  'trapz' :    trapezoidal rule
             2.  'cumtrapz' : cumulative trapezoidal
             3.  'romb' :     Romberg integration
             4.  'simps' :    Simpson's rule
             
    **Returns:**
    
       * A : integrated matrix solution.  Member of the operator class.
    """
    
    # Check user supplied inputs
    if not ctrl.number_controls == len(hamiltonians):
        raise ValueError('Bilinear dimension mismatch.')
    
    
    # Create an array to hold numeric values of integrals.  Since
    # integration is linear, we may integrate each control function
    # independently and perform the nessisary summation at the end.
    intgrl = zeros( [ctrl.number_controls , 1] )
    t = ctrl.times.flatten()
    
    if method == 'trapz':
        
        for index in range( ctrl.number_controls ):
            y = ctrl.control[:,index].flatten()
            intgrl[index] = trapz(y,t)
            
    elif method == 'cumtrapz':
        
        for index in range( ctrl.number_controls ):
            y = ctrl.control[:,index].flatten()
            intgrl[index] = cumtrapz(y,t)
        
    elif method == 'romb':
        
        for index in range( ctrl.number_controls ):
            y = ctrl.control[:,index].flatten()
            intgrl[index] = romb(y,t)
        
    elif method == 'simps':
        
        for index in range( ctrl.number_controls ):
            y = ctrl.control[:,index].flatten()
            intgrl[index] = simps(y,t)
            
    else:
        raise ValueError('Method %s not understood.' %(method))
    
    # Multiply the Hamiltonians by the required integrals
    A = 0 * hamiltonians[0]
    for index in range( ctrl.number_controls ):
        A = A + float(intgrl[index]) * hamiltonians[index]
        
    # Return the integrated matrix
    return A
    
    
def trotter( ctrl, hamiltonians ):
    """
    Solves a bilinear control system using a Trotter formula.
    
    **Forms:**
    
        ``trotter( ctrl, hamiltonians )``
        
    **Args:**
    
        * *ctrl* :   An instance of the control class.  Contains 
          time information as well as k-many control functions.
        * *hamiltonians* :  A list or array of k-many Hamiltonians.  
          The Hamiltonians must be square matrices of the same 
          dimensionality.
          
    **Returns:**
    
        * U : Solution to bilinear control problem.
    """
    
    # Initialize propagator for the control system
    dimension = sqrt( hamiltonians[0].size )
    U = operator( eye( dimension ) )
    
    # For each discrete timestep, the controls are assumed constant.
    # Then on each timestep the propagator is exp(-i H dt).  For each
    # timestep, calculate the Hamiltonian and pulse duration.
    
    for timestep in range( len(ctrl.times) - 1 ):
        
        # Calculate pulse duration
        dt = ctrl.times[ timestep + 1 ] - ctrl.times[ timestep ]
        dt = float(dt)
        
        # Construct Hamiltonian over this interval
        c = ctrl.control[timestep,:]
        H = 0 * hamiltonians[0]
        for index in range( ctrl.number_controls ):
            H = H + c[index] * hamiltonians[index]
            
        # Compute propagator for this interval
        Ut = operator( expm( -1j * H * dt ) )
        
        # Append evolution
        U = operator( Ut * U )
        
    # Return propagator
    return U


def dyson( ctrl, hamiltonians, order = 4 ):
    """
    Solves a bilinear control system using a Dyson series.  By
    default, the series is truncated at fourth order.  The series is
    not checked for convergence.
    
    **Forms:**
    
        * ``dyson( ctrl, hamiltonians )``
        * ``dyson( ctrl, hamiltonians, order = 4 )``
        
    **Args:**
    
        * *ctrl* : An instance of the control class.  Contains time
          information as well as k-many control functions.
        * *hamiltonians* : A list or array of k-many Hamiltonians.  
          The Hamiltonians must be square matrices of the same 
          dimensionality.
                 
    **Optional keys:**
    
        * order :  The order of the Dyson series.  Must be an integer.
        
    **Returns:**
    
        * U : Solution to bilinear control problem.
    """
    return NotImplemented


def magnus( ctrl, hamiltonians, order = 4 ):
    """
    Solves a bilinear control system using a Magnus expansion.  By
    default, the series is truncated at fourth order.  The expansion
    is not checked for convergence.
    
    **Forms:**
    
        * ``magnus( ctrl, hamiltonians )``
        * ``magnus( ctrl, hamiltonians, order = 4 )``
        
    **Args:**
    
        * *ctrl* : An instance of the control class.  Contains time
          information as well as k-many control functions.
        * *hamiltonians* :  A list or array of k-many Hamiltonians.  
          The Hamiltonians must be square matrices of the same 
          dimensionality.
                 
    **Optional keys:**
    
        * order :  The order of the Magnus series.  Must be an integer.
        
              
    **Returns:**
    
        * U : Solution to bilinear control problem.
    """
    return NotImplemented


def lindblad( ctrl, hamiltonians, channels ):
    return NotImplemented
