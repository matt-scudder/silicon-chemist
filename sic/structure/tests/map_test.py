"""
This class tests whether we get the right maps for particular molecules.
When you encounter a particularly problematic map, or bugs with AAMTool,
add that case to this file.
"""

import unittest

from openbabel.pybel import readstring

from sic.structure import properties

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
        self.assertDictEqual(mapping[0], correct_mapping)

if __name__ == "__main__":
    unittest.main()
