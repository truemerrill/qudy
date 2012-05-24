#!/usr/bin/env python
#
# SK1.py
#
# A demonstration of qudy.  In this example we construct a SK1 pulse
# sequence that implements the rotation R( pi/2, 0 ).  The SK1
# sequence improves the accuracy of gates in the presence of
# systematic amplitude control errors.

from qudy import *
from qudy.quantop import *
import qudy.plot as plt
"""
try:
    import matplotlib.pyplot as plt

except ImportError:
    print "Requires matplotlib for plotting."
"""

# Create a systematic error
err = error('amplitude')

# Construct the pulse sequences
theta = pi/2
U = R( theta , 0 )                  # Ideal gate
V = M( theta , 0 , err )            # Naive pulse

phi = arccos( - theta / (4*pi) ) 
SK1 = M(2*pi, -phi, err) * M(2*pi, phi, err) * M(theta, 0, err )

# Plot the controls
SK1.control.plot()

# Create a scaling plot
plt.scaling( V, U )
plt.scaling( SK1, U )

"""
# Matrix representation for the ideal gate
U_matrix = U.solve()

# Create a set of error parameters to scan over
epsilon = 10 ** arange( 0, -6, -0.1 )
raw_data = zeros( epsilon.shape )
SK1_data = zeros( epsilon.shape )

for index in range( len(epsilon) ):
    
    # Update V to a new strength of the control error.
    V.error.error_parameters = [ epsilon[index] ]
    V.update_error()
    
    # Update SK1 to a new strength of the control error.
    SK1.error.error_parameters = [ epsilon[index] ]
    SK1.update_error()
    
    # Solve the control problem.  Calculate infidelity of gate
    # relative to the ideal operation.
    V_matrix = V.solve()
    SK1_matrix = SK1.solve()
    
    raw_data[index] = infidelity( U_matrix, V_matrix )
    SK1_data[index] = infidelity( U_matrix, SK1_matrix )
   

# Create a scaling plot from infidelity data

plt.figure()    
plt.loglog( epsilon, raw_data, \
            epsilon, SK1_data )
plt.xlabel( r'$\epsilon$' )
plt.ylabel( r'Infidelity' )
"""
