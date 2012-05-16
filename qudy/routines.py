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

__all__ = ['inner_product','projection','norm','decomp','trace_distance', \
           'fidelity','infidelity','gram_schmidt','commutator',    \
           'product_operator','generate_algebra',   \
           'structure_constants','euler_decomposition']

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


def projection( U,V ):
    """
    Project an object V onto a basis vector U.
    """
    if inner_product(U,V) == 0:
        return 0 * U
    
    else:
        return inner_product(V,U) / inner_product(U,U) * U


def norm(A):
    """
    Hilbert-Schmidt norm between two matrices.
    
    **Forms:**
    
       * ``norm(A)``
       
    **Args:**
    
       * *A* : a N x N dimensional matrix
       
    **Returns:**
    
       * nrm : norm, defined as :math:`\\sqrt{ \\langle A, A \\rangle }`.
    """
    nrm = sqrt( inner_product(A,A) )
    return nrm


def decomp(A):
    """
    Decomposes a matrix U in SU(2) into coefficents (ax,ay,az) such
    that U = exp( -i/2 * (ax * X + ay * Y + az * Z) ).  This function
    does not use matrix logrithms.

    **Forms:**

       * ``decomp(A)``

    **Args:**

       * *A* : a 2 x 2 dimensional special unitary matrix

    **Returns:**
    
       * [ax,ay,az] : a list of coefficients
    """
    I = operator("1,0;0,1")
    
    # Check whether input is unitary
    #if not (A*A.H) == I:
    #    raise ValueError("Input matrix is not unitary.")

    # Check whether input is in SU(2)
    if not real( det(A) ) == 1:
        raise ValueError("Input matrix is not special unitary.")

    alpha = 2*arccos( real( A[0,0] ) )
    beta  = - sin( alpha / 2.0 )

    # Find components
    ax = alpha/beta * imag( A[0,1] )
    ay = alpha/beta * real( A[0,1] )
    az = alpha/beta * imag( A[0,0] )

    return [ax,ay,az]
    

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
    defined as :math:`F = \\min_{\\psi} \\sqrt{ \\langle \\psi | A^\\dagger B |\\psi \\rangle \\langle \\psi | B^\\dagger A | \\psi \\rangle }`.
       
    **Forms:**
    
       * ``fidelity( A, B )``
       
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
    return float( real(fidlty) )


def infidelity(A, B):
    """
    Calculates the infidelity beween to unitaries.  Infidelity is
    defined as :math`1 - F(A,B)`, where :math:`F(A,B)` is the
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


def gram_schmidt( vectors ):
    """
    From a list of vectors, generate an orthanormal basis.  If one
    element is linearly dependent on the others, then it is removed.
    """
    
    # Initialize the Gram Schmidt procedure
    e = vectors[0] / norm(vectors[0])
    basis = [ e ]
    
    # Perform Gram Schmidt procedure
    for v in vectors[1:]:

        u = v
        for e in basis:
            u = u - projection( e, v )
        
        # Normalize basis element, assuming that the norm is not equal
        # to zero.
            
        if not norm(u) == 0:
            e = u / norm(u)

            # Check to see if e is already in basis
            duplicate = False
            for b in basis:
                if not inner_product( e , b ) == 0:
                    duplicate = True

            if not duplicate:
                basis.append(e)
            
    return basis


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


def generate_algebra( hamiltonians, max_depth = 10 ):
    """
    Produces a set of skew-symmeterized Hamiltonians which completely
    span the dynamical algebra from an initial seed set.  The
    Hamiltonians act as generators, i.e. the algebra is produced by
    repeated brackets [-iH,-iH].
    """
    
    # Take the list of Hamiltonians and convert them to generators.
    generators = []
    for h in hamiltonians:
        generators.append( -1j * h )
        
    # Orthanormalize the set
    generators = gram_schmidt( generators )
    brackets = generators[:]
    
    for depth in range(1, max_depth + 1):
        
        # Measure dimension of original set
        dimension = len(brackets)
        
        # Calculate new lie brackets
        new_brackets = []

        for b in brackets:
            for g in generators:
                
                new_brackets.append( commutator(b,g) )
                
        # Orthogonalize the list of brackets to build a larger basis.
        # This throws out repeat brackets.
        brackets = gram_schmidt( brackets + new_brackets )
        
        if len(brackets) == dimension:
            
            # print "Converged at depth k = %i." %(depth)
            # print "Algebra dimension d = %i." %(dimension)
            return brackets
        
    # Must not have converged
    raise ValueError("Failed to converge at depth k = %i.\n" %(depth) + \
                     "Algebra dimension d > %i." %( len(brackets) ))
    
    return brackets


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
    Decompose a unitary in :math:`U(2)` into a set of three Euler
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
         R_x(\\gamma) R_y(\\beta) R_x(\\alpha)`.
         
    .. todo::
       
       Check stability of algorithm for rotations in the plane.  Check
       that elements U(2) are mapped correctly to equivalent elements
       in SU(2).
    """
    
    # Check to see if U is in SU(2).  If it is not, but still
    # two-dimensional, adjust the global phase so that det(U) = 1.
    # Otherwise throw an error.
    
    if not U.shape == (2,2):
        raise ValueError("Matrix must be two dimensional")
    
    if not det(U) == 1:
        
        # Adjust global phase so that determinant is 1.
        phase = 1j/2.0 * log( det(U) )
        U = exp( 1j * phase ) * U
    
    
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
