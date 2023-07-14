"""
This class tests whether we get the right maps for particular molecules.
When you encounter a particularly problematic map, or bugs with AAMTool,
add that case to this file.
"""
import unittest
import structure.properties as properties
from pybel import readstring

class MappingTest(unittest.TestCase):
    def setUp(self):
        reactants = readstring("smi","CC(O)(C)C.Cl")
        products = readstring("smi","CC(Cl)(C)C.O")
        reactants.addh()
        products.addh()
        self.reactants = reactants
        self.products = products

    def testMapping(self):
        #mapping should be:
        #1:6, 2:2, 3:3, 4:4, 5:5, 6:1
        correct_mapping = {1:6,2:2,3:3,4:4,5:5,6:1}
        mapping = properties.get_mapping(self.reactants,self.products)
        self.assertTrue(mapping == correct_mapping) #amazingly, this works. Bless Python.

if __name__ == "__main__":
    unittest.main()
