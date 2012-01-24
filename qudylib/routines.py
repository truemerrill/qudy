# ROUTINES.PY
#
# Miscellaneous subroutines for quantum control.
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

# ******************************************************
# Distance Measures                                    *
# ******************************************************

def inner_product(A, B):
    """
    Hilbert-Schmidt inner product between two matrices.
    
    **Forms:**
    
       ``inner_product(A,B)``
       
    **Args:**
    
       * *A* : a N x N dimensional matrix
       * *B* : a N x N dimensional matrix
       
    **Returns:**
    
       * ip : inner product, defined as :math:`tr( A^\dagger B )`.
    """
    ip = trace( A.H * B )
    return ip


def norm(A):
    """
    Hilbert-Schmidt norm between two matrices.
    
    **Forms:**
    
       * ``norm(A)``
       
    **Args:**
    
       * *A* : a N x N dimensional matrix
       
    **Returns:**
    
       * nrm : norm, defined as :math:`\sqrt{ \langle A, A \rangle }`.
    """
    nrm = sqrt( inner_product(A,A) )
    return nrm


def trace_distance(A, B):
    """
    Calculates the trace distance between two matrices.  Several
    definitions for the trace distance exist, here we define the
    distance between matrices :math:`A` and :math:`B` as :math:`|| A -
    B ||`.
    
    **Forms:**
       
       * `trace_distance(A,B)`
    
    **Args:**
    
       * *A* : a N x N dimensional matrix
       * *B* : a N x N dimensional matrix
       
    **Returns:**
    
       * dist : trace distance, :math:`|| A - B ||`
    """
    dist = norm( A - B )
    return dist


def fidelity(A, B):
    """
    Calculates the fidelity between two unitaries.  The fidelity is
    defined as 
    
    .. math::
       
       F = min_\psi \sqrt{\langle \psi | A^\dagger B |\psi\rangle \langle \psi | B^\dagger A | \psi \rangle }.
       
    **Forms:**
    
       * `fidelity( A, B )`
       
    **Args:**
    
       * *A* : a N x N dimensional matrix
       * *B* : a N x N dimensional matrix
       
    **Returns:**
    
       * fidlty : fidelity measure between input matrices.
    """
    [q,Q] = eig(A.H * B)
    # print q
    f = 1j* log(q)
    f = f - min(f)
    fidlty = min( abs( cos( f/2.0 ) ) )
    return real(fidlty)


def infidelity(A, B):
    """
    Calculates the infidelity beween to unitaries.  Infidelity is
    defined as :math"`1 - F(A,B)`, where :math:`F(A,B)` is the
    fidelity between matrices :math:`A` and :math:`B`.
    
    **Forms:**
    
       * `infidelity( A, B )`
       
    **Args:**
       
       * *A* : a N x N dimensional matrix
       * *B* : a N x N dimensional matrix
       
    **Returns:**
    
       * infd : infidelity measure between input matrices.
    """
    [q,Q] = eig(A.H * B)
    f = real( 1j * log(q) )
    f = f - mean(f)
    
    # Remove over rotations
    #for index in range(len(f)):
    #    if abs( f[index] ) > pi/2.0:
    #        f[index] = f[index]%(pi/2.0)
        
    dst = 2*max( sin(f/2) )
    infd = dst**2 / 2
    if infd > 1:
        infd = 2.0 - infd
    return infd




# ******************************************************
# Lie algebra basis / Generators                       *
# ******************************************************

def commutator(A, B):
    """
    Calculates the matrix commutator :math:`[A,B] = AB-BA`.
    
    **Forms:**
    
       * `commutator( A, B )`
       
    **Args:**
    
       * *A* : a N x N dimensional matrix
       * *B* : a N x N dimensional matrix
       
    **Returns:**
    
       *  : commutator between input matrices.
    """
    return A*B - B*A


def product_operator( number_qubits ):
    """
    A function to generate product operators.  The product operator
    representation is a set of orthogonal Hamiltonians which span the
    Lie algebra :math:`su(2^n)`.  This is the dynamical Lie algebra
    for a n-qubit system.
    
    **Forms:**
    
       * `product_operator( number_qubits )`
       
    **Args:**
    
       * *number_qubits* : number of qubits for Hamiltonian interactions.
         This input specifies the dimensionality of the Lie algebra, 
         :math:`4^n - 1` and the dimensionality of the Hilbert space,
         :math:`2^n`.
         
    **Raises:**
    
       * `ValueError` : An integer number of qubits is required.
         
    **Returns:**
    
       * prod_operators : A list of Hamiltonians which span the Lie algebra.
    
    See `Y Tomita et al 2010 New J. Phys. 12 015002
    <http://iopscience.iop.org/1367-2630/12/1/015002>`_ for
    conventions.
    """
    
    # Wrapper function for Pauli operators
    def s(x):
        if x == 0:
            return operator("1,0;0,1")
        elif x == 1:
            return operator("0,1;1,0")
        elif x == 2:
            return operator("0,-1j;1j,0")
        elif x == 3:
            return operator("1,0;0,-1")
        else:
            raise Exception("Could not understand input.")
    
    # Check that number_qubits is an integer value.
    if not number_qubits % 1 == 0:
        ValueError('An integer number of qubits is required.')
        
        
    # Create a list of the product operators.  Uses convention in NJP
    # 12 015002 (2010)
    prod_operators = []
    for j in range( 1, int(4**number_qubits - 1 + 1) ):
        
        eta_j = operator("1")
        for k in range( 1, int(number_qubits + 1) ):
            index = floor( (j % (4**k)) / 4**(k-1) )
            sk = s(index)
            eta_j = eta_j | sk    # Tensor product together
            
        prod_operators.append(0.5 * eta_j)
        
    return prod_operators


def generate_algebra( hamiltonians ):
    """
    Produces a set of Hamiltonains which completely span the dynamical
    algebra from an initial seed set.  Not currently implemented.
    
    .. todo::
       
       Implement algorithm to produce a set of Hamiltonians that
       spans the dynamical Lie algebra.
    """
    return NotImplemented


def structure_constants( hamiltonians ):
    """
    Calculates structure constants for dynamical Lie algebra.  Not
    currently implemented.
    
    .. todo::
       
       Implement algorithm to measure structure constants from a set
       of operators that span an algebra.
    """
    return NotImplemented


def euler_decomposition( U ):
    """
    Decompose a unitary in :math:`SU(2)` into a set of three Euler
    angles.  We use the XYX Euler angle convention.  I think it only
    works for out of plane rotations (all I care about actually).
    
    **Forms:**
    
       * `euler_decomposition( U )`
       
    **Args:**
    
       * *U* : A 2-dimensional unitary matrix.  The function does not
         check for the dimensionality or the unitarity of the matrix.
    
    **Returns:**
    
       * `[alpha,beta,gamma]` : Euler angles for the decomposition.
         The original input matrix may be reconstructed by :math:`U =
         R_x(\gamma) R_y(\beta) R_x(\alpha)`.
         
    .. todo::
       
       Check stability of algorithm for rotations in the plane.
    """
    
    H = product_operator( 1 )

    # Components of U in Pauli basis
    a1 = real( trace( U * operator("1,0;0,1") * 1/2.0 ) )
    ax = real( 1j * trace( U * H[0] ) )
    ay = real( 1j * trace( U * H[1] ) )
    az = real( 1j * trace( U * H[2] ) )

    # Calculate Euler angles
    beta = 2 * arccos( sqrt( a1**2 + ax**2 ) )
    plus = 2 * arcsin( ax / cos(beta/2.0) )  # gamma + alpha
    minus = 2 * arccos( ay / sin(beta/2.0) ) # gamma - alpha
    gamma = (plus + minus) / 2.0
    alpha = (plus - minus) / 2.0

    return [alpha,beta,gamma]
