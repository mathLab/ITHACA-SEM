
from NekPy.LibUtilities import ShapeType
from NekPy.NekMesh import Node, Element, ElmtConfig, NodeSet, Mesh, \
                          Module, ProcessModule, ModuleType,        \
                          InputModule, OutputModule
import numpy as np
import unittest
import sys


class InheritFromInputModuleTest(InputModule):
    def __init__(self, mesh):
        super(InputModule, self).__init__(mesh)

    def Process(self):
        # Call the Module functions to create all of the edges, faces and
        # composites.
        self.ProcessVertices()
        self.ProcessEdges()
        self.ProcessFaces()
        self.ProcessElements()
        self.ProcessComposites()

class InheritFromProcessModuleTest(ProcessModule):
    def __init__(self, mesh):
        super(ProcessModule, self).__init__(mesh)

    def Process(self):
        # Call the Module functions to create all of the edges, faces and
        # composites.
        self.ProcessVertices()
        self.ProcessEdges()
        self.ProcessFaces()
        self.ProcessElements()
        self.ProcessComposites()

class InheritFromOutputModuleTest(OutputModule):
    def __init__(self, mesh):
        super(OutputModule, self).__init__(mesh)

    def Process(self):
        # Call the Module functions to create all of the edges, faces and
        # composites.
        self.ProcessVertices()
        self.ProcessEdges()
        self.ProcessFaces()
        self.ProcessElements()
        self.ProcessComposites()



class TestModule(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    def _initialize_static_values(self):
        self.coord_1x   = 0.0
        self.coord_1y   = 1.0
        self.coord_2y   = 2.0
        self.coord_2x   = 3.0
        self.nx         = 11
        self.ny         = 7
        self.comp_ID    = 2
        self.shape_type = ShapeType.Triangle
        self.x_points   = np.linspace(self.coord_1x, self.coord_2x, self.nx)
        self.y_points   = np.linspace(self.coord_1y, self.coord_2y, self.ny)
        self.expDim     = 2
        self.spaceDim   = 2
        self.nummode    = 5
        self.verbose    = True

    def _create_mesh(self):
        self.mesh = Mesh()
        self.mesh.expDim   = self.expDim
        self.mesh.spaceDim = self.spaceDim
        self.mesh.nummode  = self.nummode
        self.mesh.verbose  = self.verbose

    def _initialize_nodes(self):
        self.nodes = []
        id_cnt = 0

        for y in range(self.ny):
            tmp = []
            for x in range(self.nx):
                tmp.append(Node(id_cnt, self.x_points[x], self.y_points[y], 0.0))
                id_cnt += 1
            self.nodes.append(tmp)

    def _create_triangular_elements(self):
        for y in range(self.ny - 1):
            for x in range(self.nx - 1):
                config = ElmtConfig(ShapeType.Triangle, 1, False, False)
                self.mesh.element[2].append(
                    Element.Create(
                    config,
                    [self.nodes[y][x], self.nodes[y+1][x+1], self.nodes[y+1][x]],
                    [self.comp_ID]))
                self.mesh.element[2].append(
                    Element.Create(
                    config,
                    [self.nodes[y][x], self.nodes[y][x+1], self.nodes[y+1][x+1]],
                    [self.comp_ID]))

    def _setUpMesh(self):
        self._initialize_static_values()
        self._create_mesh()
        self._initialize_nodes()
        self._create_triangular_elements()

    def setUp(self):
        self._setUpMesh()

    def testModuleProcessRuntimeError(self):
        msg = self.getCN() + "::testModuleProcessRuntimeError: "
        mod = ProcessModule(self.mesh)
        try:
            mod.Process()
            msg += "FAIL"
            print(msg)
        except RuntimeError:
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testInheritFromInputModuleTest(self):
        msg = self.getCN() + "::testInheritFromInputModuleTest: "
        mod = InheritFromInputModuleTest(self.mesh)
        try:
            mod.Process()
            # OutputModule.Create("xml", self.mesh, outfile=sys.stdout).Process()
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testInheritFromProcessModuleTest(self):
        msg = self.getCN() + "::testInheritFromProcessModuleTest: "
        mod = InheritFromProcessModuleTest(self.mesh)
        try:
            mod.Process()
            # OutputModule.Create("xml", self.mesh, outfile=sys.stdout).Process()
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testInheritFromOutputModuleTest(self):
        msg = self.getCN() + "::testInheritFromOutputModuleTest: "
        mod = InheritFromOutputModuleTest(self.mesh)
        try:
            mod.Process()
            # OutputModule.Create("xml", self.mesh, outfile=sys.stdout).Process()
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

if __name__ == '__main__':
    unittest.main()
