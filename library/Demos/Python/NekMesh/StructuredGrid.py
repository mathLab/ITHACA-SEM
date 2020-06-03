import NekPy
import sys
from NekPy.LibUtilities import ShapeType
from NekPy.NekMesh import Node, Element, ElmtConfig, NodeSet, Mesh, Module, ModuleType
import numpy as np

#
# StructuredGrid creates a 2D structured grid of triangles or quads.
#
class StructuredGrid(Module):
    def __init__(self, mesh):
        super().__init__(mesh)

        # This for now only works in 2D.
        self.mesh.spaceDim = 2
        self.mesh.expDim = 2

        # Define some configuration options for this module. AddConfigOption
        # requires 3 arguments and has one optional argument:
        # - the name of the option;
        # - the default value for that option to be used if one is not supplied
        #   via RegisterConfig;
        # - a short description, which is useful if one calls PrintConfig;
        # - optionally, if the config option is supposed to be a boolean
        #   (i.e. on or off), specify True for the fourth parameter.
        self.AddConfigOption("nx", "2", "Number of points in x direction")
        self.AddConfigOption("ny", "2", "Number of points in y direction")
        self.AddConfigOption("lx", "0", "Lower-left x-coordinate")
        self.AddConfigOption("rx", "0", "Upper-right x-coordinate")
        self.AddConfigOption("ly", "0", "Lower-left y-coordinate")
        self.AddConfigOption("ry", "0", "Upper-right y-coordinate")
        self.AddConfigOption("compid", "0", "Composite ID")
        self.AddConfigOption("shape", "Quadrilateral", "Triangular/Quadrilateral Mesh")

    def Process(self):
        # Get the input variables from our configuration options. You can use
        # either GetFloatConfig, GetIntConfig, GetBoolConfig or GetStringConfig
        # depending on the type of data that should have been input into the
        # configuration option: in the below we know that the length is a float
        # and the numbers of nodes are ints.
        coord_1x   = self.GetFloatConfig("lx")
        coord_1y   = self.GetFloatConfig("ly")
        coord_2x   = self.GetFloatConfig("rx")
        coord_2y   = self.GetFloatConfig("ry")
        nx         = self.GetIntConfig("nx")
        ny         = self.GetIntConfig("ny")
        comp_ID    = self.GetIntConfig("compid")
        shape_type = self.GetStringConfig("shape")
        x_points   = np.linspace(coord_1x, coord_2x, nx)
        y_points   = np.linspace(coord_1y, coord_2y, ny)

        nodes = []
        id_cnt = 0

        for y in range(ny):
            tmp = []
            for x in range(nx):
                tmp.append(Node(id_cnt, x_points[x], y_points[y], 0.0))
                id_cnt += 1
            nodes.append(tmp)

        if shape_type[0].lower() == "q":
            self._create_quadrilaterals(nodes, nx, ny, comp_ID)
        elif shape_type[0].lower() == "t":
            self._create_triangles(nodes, nx, ny, comp_ID)
        else:
            raise ValueError("Unknown shape type: should be quad or tri.")

        # Call the Module functions to create all of the edges, faces and
        # composites.
        self.ProcessVertices()
        self.ProcessEdges()
        self.ProcessFaces()
        self.ProcessElements()
        self.ProcessComposites()

    def _create_quadrilaterals(self, nodes, nx, ny, comp_ID):
        for y in range(ny-1):
            for x in range(nx-1):
                config = ElmtConfig(ShapeType.Quadrilateral, 1, False, False)
                self.mesh.element[2].append(
                    Element.Create(
                        config, # Element configuration
                        [nodes[y][x], nodes[y][x+1], nodes[y+1][x+1], nodes[y+1][x]], # node list
                        [comp_ID])) # tag for composite.

    def _create_triangles(self, nodes, nx, ny, comp_ID):
        for y in range(ny-1):
            for x in range(nx-1):
                config = ElmtConfig(ShapeType.Triangle, 1, False, False)
                self.mesh.element[2].append(
                    Element.Create(
                    config,
                    [nodes[y][x], nodes[y+1][x+1], nodes[y+1][x]],
                    [comp_ID]))
                self.mesh.element[2].append(
                    Element.Create(
                    config,
                    [nodes[y][x], nodes[y][x+1], nodes[y+1][x+1]],
                    [comp_ID]))

# Register our TestInput module with the factory.
Module.Register(ModuleType.Input, "StructuredGrid", StructuredGrid)

if __name__ == '__main__':
    if len(sys.argv) != 10:
        print("Usage: StructuredGrid.py nx ny lx ly rx ry compid shape outputfile.xml")
        print("")
        print("  e.g. for a 5x6 quadrilateral grid from (0,1) -> (2,3) with domain composite ID 2:")
        print("  StructuredGrid.py 5 6 0.0 1.0  2.0 3.0  2 Quad  output.xml")
        exit(1)

    # Create a 'pipeline' of the input and output modules.
    mesh = Mesh()
    mod = [
        Module.Create(ModuleType.Input, "StructuredGrid", mesh,
                      nx = sys.argv[1],
                      ny = sys.argv[2],
                      lx = sys.argv[3],
                      ly = sys.argv[4],
                      rx = sys.argv[5],
                      ry = sys.argv[6],
                      compid = sys.argv[7],
                      shape = sys.argv[8]),
        Module.Create(ModuleType.Output, "xml", mesh, outfile=sys.argv[9])
    ]

    # Print out config options that we registered in StructuredGrid, just for
    # fun and to test they registered correctly!
    print("Available options from structured grid generator:")
    mod[0].PrintConfig()

    # Process the pipeline.
    for m in mod:
        m.Process()
