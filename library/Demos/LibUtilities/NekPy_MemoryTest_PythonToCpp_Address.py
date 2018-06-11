from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList2D
import ctypes
import gc
import sys

def main():
	session_name = ["memory-test-python-to-c-address.py", "newsquare_2x2.xml"]

	session = SessionReader.CreateInstance(session_name)
	graph = MeshGraph.Read(session)
	exp = ExpList2D(session, graph)

	coords = exp.GetCoords()
	coords = coords[0]

	coords_data_address = coords.ctypes.data
	print "The memory address of data is: {}\n".format(hex(coords_data_address))

	exp.SetPhysArray(coords)
	print exp.GetPhysAddress()

if __name__ == '__main__':
    main()
