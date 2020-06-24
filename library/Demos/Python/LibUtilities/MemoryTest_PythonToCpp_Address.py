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

    coords_data_address = str(hex(coords.ctypes.data)).rstrip('L')

    print( "The memory address of data is: {}\n".format(coords_data_address))

    exp.SetPhysArray(coords)
    phys_data_address = exp.GetPhysAddress()
    print( "The memory address of array data is: {}\n".format(phys_data_address))

    if coords_data_address == phys_data_address:
        print("Test successful!")
    else:
        print("Test unsuccessful")

if __name__ == '__main__':
    main()
