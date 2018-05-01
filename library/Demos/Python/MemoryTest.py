from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList
import sys

session = SessionReader.CreateInstance(sys.argv)
graph = MeshGraph.Read(session)
exp = ExpList2D(session, graph)

print("Loaded session: %s" % session.GetSessionName())
print("Loaded MeshGraph of dimension: %d" % graph.GetMeshDimension())
