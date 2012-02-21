#!/usr/bin/env python
#
# arrays.py
#
# A demonstration of qudy.  In this example we show how Numpy
# arrays may be used to construct control instances.  
# 
# Copyright (C) 2011, 2012 True Merrill
# Georgia Institute of Technology
# School of Chemistry and Biochemistry

from qudy import *
from qudy.quantop import *

# Set up some control functions
ux = lambda t: cos( pi * t )
uy = lambda t: sin( pi * t )
uz = lambda t: 0

# Set up the time interval
dt = 1.0 / 150.0
t = arange( 0, 1 + dt, dt )

# Create the control instance
ctrl = control( ux, uy, uz, t )
ctrl.plot()

# Make a propagator.  Since this propagator is in SU(2), qudy is smart
# enough to provide the Hamiltonians we need.
U = propagator( ctrl )
print U.solve()
