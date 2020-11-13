from NekPy.NekMesh import ElmtConfig, Node, Element
from NekPy.LibUtilities import ShapeType, PointsType

import unittest

class TestElmtConfig(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    def testElmtConfigConstructor(self):
        msg = self.getCN() + "::testElmtConfigConstructor: "
        try:
            self.shapeType = ShapeType.Triangle
            self.order     = 2
            self.faceNodes = True
            self.volNodes  = True
            self.reorient  = 2
            self.edgeNodeType = PointsType.GaussLobattoLegendre
            self.faceNodeType = PointsType.GaussGaussLegendre
            self.config = ElmtConfig(self.shapeType, self.order,
                                     self.faceNodes, self.volNodes,
                                     self.reorient,  self.edgeNodeType,
                                     self.faceNodeType )
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNoDefaultConstructor(self):
        msg = self.getCN() + "::testNoDefaultConstructor: "
        try:
            default_config = ElmtConfig()
            msg += "FAIL"
        except :
            msg += "PASS"
        print(msg)


class TestElement(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    def setUp(self):
        self.x  = float(1.1)
        self.y  = float(2.2)
        self.z  = float(3.3)
        self.node_a = Node(1, self.x      , self.y      , self.z)
        self.node_b = Node(2, self.x + 1.0, self.y      , self.z)
        self.node_c = Node(3, self.x      , self.y + 1.0, self.z)
        self.config = ElmtConfig(ShapeType.Triangle, 1, False, False)
        self.comp_ID = 2
        self.element = Element.Create(self.config,
                                  [self.node_a, self.node_b, self.node_c],
                                  [self.comp_ID])

    def testElementGetId(self):
        msg = self.getCN() + "::testElementGetId: "
        try:
            # How to test GetId ??
            print("id:  ")
            print(self.element.GetId())
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testElementGetDim(self):
        msg = self.getCN() + "::testElementGetDim: "
        try:
            self.assertEqual( self.element.GetDim(), 2)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testElementGetShapeType(self):
        msg = self.getCN() + "::testElementGetShapeType: "
        try:
            self.assertEqual( self.element.GetShapeType(), ShapeType.Triangle)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testElementGetTag(self):
        msg = self.getCN() + "::testElementGetTag: "
        try:
            self.assertEqual( self.element.GetTag(), "T")
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise



if __name__ == '__main__':
    unittest.main()
