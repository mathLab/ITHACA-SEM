from NekPy.NekMesh import Mesh, Node, ElmtConfig, Element, NodeSet
from NekPy.LibUtilities import ShapeType
import unittest
import numpy as np

class TestMesh(unittest.TestCase):

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
        self.expDim     = 3
        self.spaceDim   = 3
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

    def setUp(self):
        self._initialize_static_values()
        self._create_mesh()
        self._initialize_nodes()
        self._create_triangular_elements()

    def testMeshConstructor(self):
        msg = self.getCN() + "::testMeshConstructor: "
        try:
            self.assertEqual(self.mesh.expDim   , self.expDim)
            self.assertEqual(self.mesh.spaceDim , self.spaceDim)
            self.assertEqual(self.mesh.nummode  , self.nummode)
            self.assertEqual(self.mesh.verbose  , self.verbose)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testMeshFieldAccess(self):
        msg = self.getCN() + "::testMeshFieldAccess: "
        self.nodeset = NodeSet()
        self.nodeset_def_len = 10
        for i in range(self.nodeset_def_len):
            id = i
            x  = float(i)
            y  = float(self.nodeset_def_len + i)
            z  = float(self.nodeset_def_len**2 + i)
            n  = Node(id, x, y, z)
            self.nodeset.add(n)
        try:
            self.expDim   = 2
            self.spaceDim = 2
            self.nummode  = 7
            self.verbose  = False
            self.mesh.node     = self.nodeset
            self.mesh.expDim   = self.expDim
            self.mesh.spaceDim = self.spaceDim
            self.mesh.nummode  = self.nummode
            self.mesh.verbose  = self.verbose
            for node in self.nodeset:
                self.assertTrue(node in self.mesh.node)
            self.assertEqual(self.mesh.expDim   , self.expDim)
            self.assertEqual(self.mesh.spaceDim , self.spaceDim)
            self.assertEqual(self.mesh.nummode  , self.nummode)
            self.assertEqual(self.mesh.verbose  , self.verbose)
            try:
                self.mesh.element = None
            except AttributeError:
                msg += "PASS"
                print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

if __name__ == '__main__':
    unittest.main()
