import sys
from NekPy.NekMesh import Mesh, ProcessModule, OutputModule

mesh          = Mesh()
mesh.expDim   = 3
mesh.spaceDim = 3
mesh.nummode  = 5
mesh.verbose  = True

if len(sys.argv) != 3:
    print('Usage: LoadCAD.py input.stp output.xml')
    print('')
    print('Generates a surface mesh of a 3D input STEP file using fixed parameters:')
    print('  min_delta = 0.04, max_delta=0.2, epsilon=0.02')

# Load the CAD file
ProcessModule.Create("loadcad", mesh, filename=sys.argv[1], verbose=True).Process()

# Load the octree
ProcessModule.Create("loadoctree", mesh, mindel=0.04, maxdel=0.2, eps=0.02).Process()

# Create a surface mesh
ProcessModule.Create("surfacemesh", mesh).Process()

# Output a 2D manifold mesh
mesh.expDim = 2

# Create a high-order surface
ProcessModule.Create("hosurface", mesh).Process()

# Dump out elemental Jacobians
ProcessModule.Create("jac", mesh, list=True).Process()

# Dump out the surface mesh.
OutputModule.Create("xml", mesh, sys.argv[2], test=True).Process()
