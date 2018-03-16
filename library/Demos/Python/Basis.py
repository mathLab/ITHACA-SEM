##
## NekPy example: basic interactions with the Points and Basis classes.
##
## This example shows some of the functionality of the Points and Basis
## wrappers. It constructs the Modified C^0 basis (found in the Karniadakis &
## Sherwin book) on a set of Gauss-Lobatto-Legendre points on the standard
## segment [-1,1]. It then constructs a mass matrix and uses this to perform a
## Galerkin L^2 projection of the function u(x) = sin(x).
##
## Functionally this is equivalent to the StdProject demo in Nektar++. You can
## check the output is correct by running the Python script StdProject.py, which
## uses the Nektar++ StdRegions classes directly.
##

# Import necessary LibUtilities components
from math import sin
from NekPy.LibUtilities import PointsKey, Points, Basis, BasisKey, PointsType, BasisType
import numpy as np

# Create a PointsKey and a Points object for numPts GLL points
numPts = 10
ptsKey = PointsKey(numPts, PointsType.GaussLobattoLegendre)
pts    = Points.Create(ptsKey)

# Print out some basic information
print("Created PointsKey with %d points" % pts.GetNumPoints())
print("Points locations:")
print(pts.GetZ())

# Test our integration points on something nice like sin(x), which has exact
# integral zero over [-1, 1]
z, w = pts.GetZW()
integral = sum([ sin(z[i]) * w[i] for i in range(0, numPts) ])
print("\nIntegral of sin(x) over [-1,1] = %g", integral)

# Create the modified C^0 basis on our integration points
numModes = numPts-1
basisKey = BasisKey(BasisType.Modified_A, numModes, ptsKey)
basis    = Basis.Create(basisKey)

print("\nCreated Basis with %d modes" % basis.GetNumModes())

# Get the basis evaluated at the points. We take a copy since GetBdata returns a
# read-only array. Note also the use of "order='F'" -- this is to store the
# resulting matrix in column-major format to align with Nektar++ storage.
bData = basis.GetBdata().copy().reshape(numPts, numModes, order='F')

# Form a mass matrix = \int \phi_i \phi_j dx
massMat = np.zeros((numModes, numModes))
for i in range(0, numModes):
    for j in range(0, numModes):
        massMat[i,j] = sum([ bData[k,i] * bData[k,j] * w[k] for k in range(0, numPts) ])

# Project sin(x)
rhs = np.zeros(numModes)
for i in range(0, numModes):
    rhs[i] = sum([sin(z[j]) * bData[j,i] * w[j] for j in range(0, numPts) ])

print("\nCalculated projection coefficients:")
print(np.linalg.solve(massMat, rhs))
