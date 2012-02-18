# QUANTOP.PY
# 
# A matrix-like class to perform quantum mechanics
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


# Import required libraries.  Make simple 'shortcut' copies of routines
# used in linear algebra.  This makes each of these routines superficially
# part of the QUANTOP python module.
from numpy import Inf, NaN, isinf, isnan, isreal, isscalar, \
     ones, eye, pi, sin, cos, tan, arcsin, arccos, arctan, exp, zeros, \
     trace, matrix, array, real, imag, conj, arange, floor, ceil, mean, \
     resize, diag, set_printoptions, hstack, vstack, copy, flipud
from numpy.linalg import eig, det
from scipy.linalg import kron, logm, norm, inv, sqrtm, eigh
from scipy.linalg import expm2 as expm
from scipy import sqrt, log, log10
from warnings import warn

# Set print options to make print formatting pretty.  The numerical
# precision of the calculation is not affected.
set_printoptions(precision=4,suppress=True,linewidth=80)

# Functions to change the printing options on-the-fly
def format_long():
    """
    Function to interactivly change print options within a python
    shell.  Inhibits suppression of long decimals.
    """
    set_printoptions(suppress = False)


def format_short():
    """
    Function to interactivly change print options within a python
    shell.  Sets suppression of long decimals.
    """
    set_printoptions(suppress = True)


class operator(matrix):
    """
    class for quantum operators.
    
    The operator class inherits all functionality of the numpy.matix
    class.  Inputs may be arrays or strings.
    
    **Forms:**
    
       *Arrays*
          
          1.  ``operator( arr )``
          2.  ``operator( array( [[0,1],[1,0]] ) )``
          
       *Strings*
       
          1.  ``operator( string )``
          2.  ``operator( "0,1;1,0" )``
          
    **Args:**
    
       * *arr* : An array object
       * *string* : A string object
    
    It is well known that most quantum operators may be represented as
    matrices in an appropriate basis.  The operator class convinently
    includes several new methods relevant to quantum control problems.
    
    **Methods:**
          
       1. ``__or__`` : the bitwise or function has been redefined to
          now mean a Kronecker product.  For two operator objects ``x``
          and ``y``, now ``x|y`` is their Kronecker product.
       
       2. ``__call__( t )`` : Operator objects are callable, like a 
          function.  By default, __call__ returns self, but this 
          allows complicated matrix functions to be built using a
          lambda function.  
    """

    def __call__(self,t):
        # Make operators callable, like a function
        return self
    
    def kron(self,x,y):
        """
         Matrix Kroneker product between objects ``x`` and ``y``.
         Symbolically this is expressed as :math:`x \\otimes y`.
         
         Also, the bitwise or function has been redefined to now mean
         a Kronecker product.  For two operator objects ``x`` and
         ``y``, now ``x|y`` is their Kronecker product.
         """
        prod = operator( kron(x,y) )
        return prod

    def __or__(self,y):
        # Overwrite the standard bit-wise or and instead implement
        # the tensor product.
        return self.kron(self,y)

    def __repr__(self):
        s = repr(self.__array__()).replace('array', 'operator')
        # now, 'operator' has 8 letters, and 'array' 5, so the columns
        # don't line up anymore. We need to add three spaces.
        
        l = s.splitlines()
        for i in range(1, len(l)):
            if l[i]:
                l[i] = '   ' + l[i]
        return '\n'.join(l)

    def __eq__(self, x):
        # Compare each element of two arrays.  If each element
        # matches, then return True,  else return false.
        if isinstance(x, operator):
            boolean = ((self - x) == 0).all()
            return boolean
        
        else:
            return matrix.__eq__(self,x)
            
    def tomatlab(self, precision = 8):
        """
        Return a string representation of the operator, compatible
        with Matlab.
        """
        (rows, cols) = self.shape
        s = ''

        s = s + '[ ... \n'
        for r in range(rows):
            for c in range(cols):

                s2 = '\t% .'+str(precision)+ \
                     'e + % .'+str(precision)+ \
                     'e * 1j'

                s2 = s2 %( real(self[r,c]) , imag(self[r,c]) )
                s = s + s2

                if not c == cols-1:
                    s = s + ','

            s = s + '; ... \n'
        
        s = s + ']'
        return s

    def totext(self, filename, precision = 8):
        """
        save a string representation of the operator to a text file
        """
        contents = self.tomatlab(precision)
        f = open(filename,"w")
        f.writelines(contents)
        f.close()
