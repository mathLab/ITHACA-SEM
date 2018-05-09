from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList2D
import sys

session_name = ["MemoryTest.py", "newsquare_2x2.xml"]

session = SessionReader.CreateInstance(session_name)
graph = MeshGraph.Read(session)
exp = ExpList2D(session, graph)

print("Loaded session: %s" % session.GetSessionName())
print("Loaded MeshGraph of dimension: %d" % graph.GetMeshDimension())

coords = exp.GetCoords()
print coords[0]

exp.SetPhysArray(coords[0])