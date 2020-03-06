#
# CurvedGeom.py: Test manual construction of a single deformed quadrilateral
# element.
#

from NekPy.LibUtilities import PointsKey, BasisKey, PointsType, BasisType
from NekPy.StdRegions import MatrixType, ConstFactorMap, ConstFactorType
from NekPy.SpatialDomains import PointGeom, SegGeom, QuadGeom, Curve
from NekPy.LocalRegions import QuadExp, MatrixKey
import numpy as np

# Construct points and segments for the unit square [0,1]^2
p = [
    PointGeom(2, i, pt[0], pt[1], 0.0)
    for i, pt in enumerate([ (0,0), (1,0), (1,1), (0,1) ])
]
s = [ SegGeom(i, 2, [p[i], p[(i+1) % 4]]) for i in range(0,4) ]

# Create a quadratic curve along the base edge.
curve = Curve(1, PointsType.PolyEvenlySpaced)
curve.points = [ p[0], PointGeom(2, 4, 0.5, -0.2, 0.0), p[1] ]
s[0] = SegGeom(i, 2, [p[0], p[1]], curve)

# Construct the geometry object for the square.
geom = QuadGeom(0, s)
ptsKey = PointsKey(9, PointsType.GaussLobattoLegendre)
basisKey = BasisKey(BasisType.Modified_A, 8, ptsKey)

# Construct a LocalRegions representation.
quadExp = QuadExp(basisKey, basisKey, geom)

# Create matrix key for Helmholtz matrix.
factors = ConstFactorMap()
factors[ConstFactorType.FactorLambda] = 5.0

# Construct a Helmholtz matrix.
helmMatKey = MatrixKey(MatrixType.Helmholtz, quadExp.DetShapeType(), quadExp, factors)
helmMat = quadExp.GenMatrix(helmMatKey)

# Construct Laplacian and mass matrices.
lapMatKey  = MatrixKey(MatrixType.Laplacian, quadExp.DetShapeType(), quadExp)
massMatKey = MatrixKey(MatrixType.Mass,      quadExp.DetShapeType(), quadExp)
lapMat     = quadExp.GenMatrix(lapMatKey)
massMat    = quadExp.GenMatrix(massMatKey)

# Test the L^inf error between H and a manually-constructed L+lambda*M
# equivalent.
print("L infinity error: %.6e" % np.abs(helmMat - (lapMat + 5.0 * massMat)).sum())
