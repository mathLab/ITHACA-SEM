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
    session_name = ["memory-test-python-to-c-address.py", "newsquare_2x2.xml"]
    expected_test_outcome = [1, 2, 1, -9.187094450483263e-15]
    actual_test_outcome = []

    session = SessionReader.CreateInstance(session_name)
    graph = MeshGraph.Read(session)
    exp = ExpList2D(session, graph)

    print("Loaded session: %s" % session.GetSessionName())
    print("Loaded MeshGraph of dimension: %d\n" % graph.GetMeshDimension())

    print("Retrieving coordinates...")
    coords = exp.GetCoords()
    coords_subs = coords[1]
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

    print("Substituting the coordinates in m_phys...")
    exp.SetPhysArray(coords_subs)
    print("Substitution completed.")
    print("Reference count for expansion coordinates: %d\n" % get_refcount(coords_address))
    actual_test_outcome.append(get_refcount(coords_address))

    print("Attempting to access the coordinates in Python...")
    print("The sum of coordinates array entires is: %r" % sum(coords))
    actual_test_outcome.append(sum(coords))

    if np.allclose(actual_test_outcome, expected_test_outcome):
        print("Test successful!")
    else:
        print("Test unsuccessful")

if __name__ == '__main__':
    main()
