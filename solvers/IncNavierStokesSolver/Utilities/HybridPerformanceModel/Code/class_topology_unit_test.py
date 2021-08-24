# Andrew Gloster
# Summer 2016
# Python Version of Communication Model Implementation
# Unit Testing of Topology Class

#------------------------------------
# Import relevant modules
#------------------------------------

import unittest

#------------------------------------
# Import relevant classes
#------------------------------------

from class_topology import Topology

#------------------------------------
# Test Cases
#------------------------------------

class Testing(unittest.TestCase):

# Check we get the correct total processes
	def test0(self):
		test_case = Topology(5, 4, 10, 2)
		self.assertEqual(test_case.PROC_TOT, 20)

# Check for correct basic core labeling
	def test1(self):
		test_case = Topology(2, 5, 10, 2)
		self.assertEqual(test_case.Core, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

# Check for basic socket labeling
	def test2(self):
		test_case = Topology(2, 5, 5, 1)
		self.assertEqual(test_case.Socket, [0, 0, 0, 0, 0, 1, 1, 1, 1, 1])

# Check for basic node labeling
	def test3(self):
		test_case = Topology(2, 5, 2, 1)
		self.assertEqual(test_case.Node, [0, 0, 1, 1, 2, 2, 3, 3, 4, 4])

# Check neighbours that are identical
	def test4(self):
		test_case = Topology(2, 5, 5, 1)
		self.assertEqual(test_case.Check_Neighbour(0, 0), 0)

# Check neighbours that are on same socket
	def test5(self):
		test_case = Topology(2, 5, 5, 2)
		self.assertEqual(test_case.Check_Neighbour(0, 1), 1)	

# Check neighbours that are on same node
	def test6(self):
		test_case = Topology(2, 5, 5, 2)
		self.assertEqual(test_case.Check_Neighbour(0, 8), 2)

# Check neighbours that are on different nodes
	def test8(self):
		test_case = Topology(4, 5, 2, 2)
		self.assertEqual(test_case.Check_Neighbour(0, 9), 3)


