from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList2D
import ctypes
import gc
import sys
import numpy as np

def get_refcount(coords_address):
    gc.collect()
    return ctypes.c_long.from_address(coords_address).value

def main():
    session_name = ["NekPy_ReferenceTest_PythonDeleteFirst.py", "newsquare_2x2.xml"]
    expected_test_outcome = [1, 2, 1, 3.469446951953614e-17]
    actual_test_outcome = []

    session = SessionReader.CreateInstance(session_name)
    graph = MeshGraph.Read(session)
    exp = ExpList2D(session, graph)

    print("Loaded session: %s" % session.GetSessionName())
    print("Loaded MeshGraph of dimension: %d\n" % graph.GetMeshDimension())

    print("Retrieving coordinates...")
    coords = exp.GetCoords()
    coords = coords[0]
    coords_address = id(coords)

    print("Retrieved coordinates.")
    print("Reference count for expansion coordinates: %d\n" % get_refcount(coords_address))
    actual_test_outcome.append(get_refcount(coords_address))

    print("Setting PhysArray (exp.SetPhysArray())...")
    exp.SetPhysArray(coords)
    print("exp.SetPhysArray() completed.")
    print("Reference count for expansion coordinates: %d\n" % get_refcount(coords_address))
    actual_test_outcome.append(get_refcount(coords_address))

    print("Deleting coordinates in Python...")
    del coords

    print("Coordinates deleted in Python.")
    print("Reference count for expansion coordinates: %d\n" % get_refcount(coords_address))
    actual_test_outcome.append(get_refcount(coords_address))

    # Ensure we are in PhysState to avoid ASSERT inside ExpList.
    exp.SetPhysState(True)

    print("Attempting to calculate the integral...")
    print("Integral calculated to be: %r" % exp.PhysIntegral())
    actual_test_outcome.append(exp.PhysIntegral())

    if np.allclose(actual_test_outcome, expected_test_outcome):
        print("Test successful!")
    else:
        print("Test unsuccessful")

if __name__ == '__main__':
    main()
