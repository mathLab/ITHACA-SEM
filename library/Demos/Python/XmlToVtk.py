from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList2D
import sys

if len(sys.argv) < 3:
    print("Usage: python XmlToVtk.py input1.xml input2.xml ... output.vtu")
    exit(1)

# Load up session and create ExpList2D
session = SessionReader.CreateInstance(sys.argv[:-1])
graph = MeshGraph.Read(session)
exp = ExpList2D(session, graph)

print("Loaded %s with %d elements" % (session.GetSessionName(), exp.GetExpSize()))

# Write a VTK file
exp.WriteVTK(sys.argv[-1])

print ("Written: %s" % sys.argv[-1])
