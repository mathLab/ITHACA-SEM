##
## NekPy example: basic interactions with the StdRegions classes
##
## This example constructs some simple StdRegions shapes and shows some of the
## operations that can be performed on the wrapper classes.
##

# Import necessary LibUtilities and StdRegions components
from NekPy.LibUtilities import PointsKey, Points, Basis, BasisKey, PointsType, BasisType
from NekPy.StdRegions import StdSegExp, StdQuadExp

# Other Python imports
from math import sin, cos
import numpy as np

# Create a PointsKey and a Points object for numPts GLL points. Create a
# Modified_A basis for the segment.
numPts   = 10
numModes = 9
ptsKey   = PointsKey(numPts, PointsType.GaussLobattoLegendre)
basisKey = BasisKey(BasisType.Modified_A, numModes, ptsKey)

# Create StdSegExp
seg = StdSegExp(basisKey)

# Use GetCoords to get coordinates of the segment. Note GetCoords always returns
# a tuple.
func = np.sin(seg.GetCoords()[0])

# Project sin(x). You can check the output here with the manual implementation
# in Basis.py.
print("Coefficients of projection of sin(x):")
print(seg.FwdTrans(func))

# Now let's create a quad on a tensor product of the points and basis above.
quad = StdQuadExp(basisKey, basisKey)

# Calculate the integral of the function cos(x)cos(y) on [-1,1]^2
xi, yi   = quad.GetCoords()
func     = np.cos(xi) * np.cos(yi)
integral = quad.Integral(func)

# Print details of our calculation.
print("Integral of cos(x) * cos(y) over the square [-1,1]^2")
print("   calculated = %g" % integral)
print("   error      = %g" % abs(4.0 * sin(1) * sin(1) - integral))
