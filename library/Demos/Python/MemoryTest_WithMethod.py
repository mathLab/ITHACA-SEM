from NekPy.LibUtilities import SessionReader
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList2D
import sys
import ctypes
import gc

def get_refcount(coords_address):
	gc.collect()
	# return len(gc.get_referrers(ctypes.cast(coords_address, ctypes.py_object).value))
	return ctypes.c_long.from_address(coords_address).value

# def namestr(obj, namespace):
#     return [name for name in namespace if namespace[name] is obj]

def test(coords, coords_address, exp):
	print("Jumped into test() method...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address))

	exp.SetPhys(coords)
	print("exp.SetPhys() completed...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address))

	del coords
	print("del coords completed...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address))


def main(): 
	session_name = ["MemoryTest.py", "newsquare_2x2.xml"]

	session = SessionReader.CreateInstance(session_name)
	graph = MeshGraph.Read(session)
	exp = ExpList2D(session, graph)

	print("Loaded session: %s" % session.GetSessionName())
	print("Loaded MeshGraph of dimension: %d" % graph.GetMeshDimension())

	coords = exp.GetCoords()
	coords = coords[0]

	# print coords

	coords_address = id(coords)
	print("Got coordinates...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address))

	test(coords, coords_address, exp)
	
	print("Jumped out of test() method...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address))
	
	# print exp.GetPhys()
	print exp.PhysIntegral()
	del exp
	print("del exp completed...")
	print("Reference count for expansion coordinates: %d" % get_refcount(coords_address))

	
	# print ctypes.cast(coords_address, ctypes.py_object).value

if __name__ == '__main__':
    main()