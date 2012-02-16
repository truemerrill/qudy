#!/usr/bin/env python
#
# Insert header describing this example.

from qudy import *
from qudy.quantop import *

# Create a systematic error
err = error('amplitude')

# Construct the pulse sequences
theta = pi/2
U = R( theta , 0 )                  # Ideal gate
V = M( theta , 0 , err )            # Naive pulse

phi = arccos( - theta / (4*pi) ) 
BB1 = M(pi, phi, err) * M(2*pi, 3*phi, err) * \
      M(pi, phi, err) * M(theta, 0, err )

# Plot the controls
BB1.controls.plot()


# Create a set of error parameters to scan over
epsilon = 10 ** arange( 0, -6, -0.1 )
raw_data = zeros( epsilon.shape )
BB1_data = zeros( epsilon.shape )

