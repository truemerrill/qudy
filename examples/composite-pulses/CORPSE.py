#!/usr/bin/env python
#
# CORPSE.py
#
# A demonstration of qudy.  In this example we construct a CORPSE pulse
# sequence that implements the rotation R( pi/2, 0 ).  The CORPSE
# sequence improves the accuracy of gates in the presence of
# systematic amplitude control errors.


from qudy import *
from qudy.quantop import *

try:
    import matplotlib.pyplot as plt

except ImportError:
    print "Requires matplotlib for plotting."


# Create a systematic error
err = error('detuning')

# Construct the pulse sequences
theta = pi/2
U = R( theta , 0 )                  # Ideal gate
V = M( theta , 0 , err )            # Naive pulse

n1 = 1; n2 = 1; n3 = 0;
theta1 = 2*pi*n1 + theta/2.0 - arcsin( sin(theta/2.0) / 2.0 )
theta2 = 2*pi*n2 - 2*arcsin( sin(theta/2.0) / 2.0 )
theta3 = 2*pi*n3 + theta/2.0 - arcsin( sin(theta/2.0) / 2.0 )    
CORPSE = M(theta3,0,err) * M(-theta2,0,err) * M(theta1,0,err)

# Plot the controls
# CORPSE.control.plot()

# Matrix representation for the ideal gate
U_matrix = U.solve()

# Create a set of error parameters to scan over
epsilon = 10 ** arange( 0, -6, -0.1 )
raw_data = zeros( epsilon.shape )
CORPSE_data = zeros( epsilon.shape )

for index in range( len(epsilon) ):
    
    # Update V to a new strength of the control error.
    V.error.error_parameters = [ epsilon[index] ]
    V.update_error()
    
    # Update CORPSE to a new strength of the control error.
    CORPSE.error.error_parameters = [ epsilon[index] ]
    CORPSE.update_error()
    
    # Solve the control problem.  Calculate infidelity of gate
    # relative to the ideal operation.
    V_matrix = V.solve()
    CORPSE_matrix = CORPSE.solve()
    
    raw_data[index] = infidelity( U_matrix, V_matrix )
    CORPSE_data[index] = infidelity( U_matrix, CORPSE_matrix )
   

# Create a scaling plot from infidelity data
plt.loglog( epsilon, raw_data, \
            epsilon, CORPSE_data )
plt.xlabel( r'$\epsilon$' )
plt.ylabel( r'Infidelity' )
plt.show()
    
