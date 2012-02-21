#!/usr/bin/env python
#
# X.PY
#
# A demonstration of qudy.  In this example we construct a "X" gate
# for a sigle qubit.
# 
# Copyright (C) 2011, 2012 True Merrill
# Georgia Institute of Technology
# School of Chemistry and Biochemistry

import qudy

# Create the propagator from a simple rotation about the X axis.
# Specifically, this will produce the gate -iX, which differes from X
# by an unimportant global phase.
U = qudy.rotation( qudy.quantop.pi , 0 )
X = qudy.quantop.operator( "0, 1; 1, 0" )

# Print the matrix
print U.solve()

# Print gate fidelity
print "fidelity : %.3f" %( qudy.fidelity( X , U.solve() ))
