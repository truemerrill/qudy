#!/usr/bin/env python
#
# H.PY
#
# A demonstration of qudy.  In this example we construct a "H" gate
# for a sigle qubit.
# 
# Copyright (C) 2011, 2012 True Merrill
# Georgia Institute of Technology
# School of Chemistry and Biochemistry

from qudy import *
from qudy.quantop import *

# Create the propagator from rotations about the X, Y, and X axis.
# Use an euler decomposition to solve for the rotation angles.
H = operator("1, 1; 1, -1") / sqrt(2)

# Find an equivalent operation in SU(2) (adjust global phase).  This
# appears to be nessisary for the euler decomposition to work.

angles = euler_decomposition( H )
U = rotation( angles[0], 0 ) * \
    rotation( -angles[1], pi/2 ) * \
    rotation( angles[2], 0 )

# Print the matrix
print U.solve()

# Print gate fidelity
print "fidelity : %.3f" %( fidelity( H, U.solve() ))
