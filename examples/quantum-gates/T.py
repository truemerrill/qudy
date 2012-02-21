#!/usr/bin/env python
#
# T.PY
#
# A demonstration of qudy.  In this example we construct a "T" gate
# for a sigle qubit.
# 
# Copyright (C) 2011, 2012 True Merrill
# Georgia Institute of Technology
# School of Chemistry and Biochemistry

from qudy import *
from qudy.quantop import *

# Create the propagator from rotations about the X, Y, and X axis.
# Use an euler decomposition to solve for the rotation angles.  

T = operator([[1, 0],[ 0, exp( 1j * pi/4)]])

# Find an equivalent operation in SU(2) (adjust global phase).  This
# appears to be nessisary for the euler decomposition to work.
N = sqrt(T.size)
phase = 1j/N * log( det(T) )
P =  operator(expm( 1j * phase * operator(eye(N)) ))
T = P*T


angles = euler_decomposition( T )
U = rotation( angles[0], 0 ) * \
    rotation( -angles[1], pi/2 ) * \
    rotation( angles[2], 0 )

# Print the matrix
print U.solve()

# Print gate fidelity
print "fidelity : %.3f" %( fidelity( T, U.solve() ))
