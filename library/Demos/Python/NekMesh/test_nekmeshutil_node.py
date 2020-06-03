from NekPy.NekMesh import Node, NodeSet
import unittest

class TestNode(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    def setUp(self):
        self.id = int(1)
        self.x  = float(1.1)
        self.y  = float(2.2)
        self.z  = float(3.3)
        self.node = Node(self.id, self.x, self.y, self.z)

    def testNodeConstructor(self):
        msg = self.getCN() + "::testNodeConstructor: "
        try:
            self.assertEqual(self.node.id, self.id)
            self.assertEqual(self.node.x,  self.x)
            self.assertEqual(self.node.y,  self.y)
            self.assertEqual(self.node.z,  self.z)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeGetID(self):
        msg = self.__class__.__name__ + "::testNodeGetID: "
        try:
            self.assertEqual(self.node.GetID(), self.id)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeSetID(self):
        msg = self.getCN() + "::testNodeSetID: "
        self.id = 2
        self.node.SetID(self.id)
        try:
            self.assertEqual(self.node.GetID(), self.id)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeDistance(self):
        msg = self.getCN() + "::testNodeDistance: "
        node = Node(self.id + 1, self.x + 1, self.y + 2, self.z + 2)
        distance = self.node.Distance(node)
        expected_distance = 3.0
        try:
            self.assertAlmostEqual(distance, expected_distance)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeGetLoc(self):
        msg = self.getCN() + "::testNodeGetLoc: "
        expected_locations = [self.x, self.y, self.z]
        locations = self.node.GetLoc()
        try:
            for i in range(len(locations)):
                self.assertEqual(locations[i], expected_locations[i])
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeAbs2(self):
        msg = self.getCN() + "::testNodeAbs2: "
        expected_abs2 = self.x**2 + self.y**2 + self.z**2
        abs2 = self.node.abs2()
        try:
            self.assertAlmostEqual(abs2, expected_abs2)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeFieldAccess(self):
        msg = self.getCN() + "::testNodeFieldAccess: "
        id = int(2)
        x  = float(11.1)
        y  = float(12.2)
        z  = float(13.3)
        try:
            self.node.id = id
            self.node.x  = x
            self.node.y  = y
            self.node.z  = z
            self.assertEqual(self.node.id, id)
            self.assertEqual(self.node.x,  x)
            self.assertEqual(self.node.y,  y)
            self.assertEqual(self.node.z,  z)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise


class TestNodeSet(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    def setUp(self):
        self.nodeset = NodeSet()
        self.nodeset_def_len = 10
        for i in range(self.nodeset_def_len):
            id = i
            x  = float(i)
            y  = float(self.nodeset_def_len + i)
            z  = float(self.nodeset_def_len**2 + i)
            n  = Node(id, x, y, z)
            self.nodeset.add(n)

    def testNodeSet__len__(self):
        msg = self.getCN() + "::testNodeSet__len__: "
        slen = len(self.nodeset)
        try:
            self.assertEqual(slen, self.nodeset_def_len)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeSetClear(self):
        msg = self.getCN() + "::testNodeSetClear: "
        self.nodeset.clear()
        slen = len(self.nodeset)
        try:
            self.assertEqual(slen, 0)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeSet__iter__(self):
        msg = self.getCN() + "::testNodeSet__iter__: "
        try:
            for node in self.nodeset:
                self.assertTrue(node.GetID() <= self.nodeset_def_len)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeSet__contains__(self):
        msg = self.getCN() + "::testNodeSet__contains__: "
        try:
            for node in self.nodeset:
                self.assertTrue(node in self.nodeset)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

    def testNodeSetAdd(self):
        msg = self.getCN() + "::testNodeSetAdd: "
        new_node = Node(len(self.nodeset) + 1, 1.0, 2.0, 3.0)
        self.nodeset.add(new_node)
        try:
            self.assertTrue(new_node in self.nodeset)
            msg += "PASS"
            print(msg)
        except:
            msg += "FAIL"
            print(msg)
            raise

if __name__ == '__main__':
    unittest.main()
