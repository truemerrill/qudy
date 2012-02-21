#!/usr/bin/env python
#
# control-Z.PY
#
# A demonstration of qudy.  In this example we construct a "H" gate
# for two qubits.
# 
# Copyright (C) 2011, 2012 True Merrill
# Georgia Institute of Technology
# School of Chemistry and Biochemistry

from qudy import *
from qudy.quantop import *


# Construct a control-Z matrix, to act as a reference
CZ = operator( [[1,0,0,0], \
                [0,1,0,0], \
                [0,0,1,0], \
                [0,0,0,-1]] )

# Hamiltonians we will need:
[Hx,Hy,Hz] = product_operator(1)
I = operator( eye(2) )
hamiltonians = [ (Hz|I) , (I|Hz) , 2*(Hz|Hz) ] 

# Controls we will need
ARR = array([[ -1, -1, 1 , 0   ], \
             [ -1, -1, 1 , pi/2]] )
c = control(ARR)

# Construct a propagator.  This gate is equivalent to control-Z up to
# a global phase.
U = propagator( c, hamiltonians )

# Print the matrix
print U.solve()

# Print gate fidelity
print "fidelity : %.3f" %( fidelity( CZ , U.solve() ))
