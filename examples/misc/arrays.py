#!/usr/bin/env python
#
# functions.py
#
# A demonstration of qudy.  In this example we show how Python
# functions may be used to construct control instances.  
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

# Create an array from the control functions.  We may always specify
# control functions in terms of a control array.
ARR = hstack(( array([ map(ux,t) ]).transpose() , \
               array([ map(uy,t) ]).transpose() , \
               array([ map(uz,t) ]).transpose() , \
               array([t]).transpose() ))

# Create the control instance
ctrl = control(ARR)

