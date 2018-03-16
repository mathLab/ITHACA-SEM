from NekPy.LibUtilities import SessionReader, BasisKey, PointsKey, PointsType, BasisType
from NekPy.SpatialDomains import MeshGraph, QuadGeom
import NekPy.StdRegions
from NekPy.LocalRegions import QuadExp
import sys
import numpy as np

session = SessionReader.CreateInstance(sys.argv)
graph = MeshGraph.Read(session)

print("Loaded session: %s" % session.GetSessionName())
print("Loaded MeshGraph of dimension: %d" % graph.GetMeshDimension())

# Get a map of all quadrilateral geometries
quad_geoms = graph.GetAllQuadGeoms()

print("Mesh contains %d quads" % len(quad_geoms))

for quaditer in quad_geoms:
    i = quaditer.key()
    quad = quaditer.data()

    # Print some basic information about each quad: it's ID, coordinate
    # dimension and the number of modes in its geometric representation
    # (i.e. NUMPOINTS in the CURVED tags).
    sys.stdout.write("\n%15s : %d\n" % ("ID", quad.GetGlobalID()))
    sys.stdout.write("%15s : %d\n" % ("Coord. Dim", quad.GetCoordim()))
    sys.stdout.write("%15s : %d, %d\n" % (
        "# modes", quad.GetBasis(0).GetNumModes(), quad.GetBasis(1).GetNumModes()))

    # Fill the geometry with its curvature and generate the geometric factors
    quad.FillGeom()
    quad.GenGeomFactors()

    # Determine if this quad contains the point (0.5, 0.5)?
    sys.stdout.write("%15s : %s\n" % (
        "Contains point", "yes" if quad.ContainsPoint(np.array([0.25, 0.25])) else "no"))

    # Create a LocalRegions expansion of this quad and print out the box that
    # bounds it.
    basisKey = BasisKey(BasisType.Modified_A, 7,
                        PointsKey(8, PointsType.GaussLobattoLegendre))
    quadExp = QuadExp(basisKey, basisKey, quad)
    x, y    = quadExp.GetCoords()
    sys.stdout.write("%15s : [%f,%f] -> [%f,%f]\n" %
                     ("Bounding box", x.min(), y.min(), x.max(), y.max()))
