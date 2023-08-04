"""
Tests whether C-L can be identified.
"""

import unittest

from openbabel.pybel import readstring,Smarts

from sic.segmentation import sinks

class CLIdentificationTest(unittest.TestCase):
    def setUp(self):
        self.cl = ["CCl","CC(C)O","CC(C)(C)O","COC(C)C","CCC(Cl)Br"]
        self.not_cl = ["C","CC(C)","CC(C)(C)","C=O","O=CCl"]
        self.SINKS = sinks.SINKS

    def testCLMatch(self):
        """
        Tests whether C-L are identified correctly as sinks when they are,
        and NOT when there's no C-L.
        """

        smarts = Smarts(self.SINKS["C-L"])
        for match in self.cl[0:3]:
            mol = readstring("smi",match)
            mol.addh()
            results = smarts.findall(mol)
            self.assertEquals(len(results), 1) #only one set of atoms 
            self.assertEquals(len(results[0]), 2) #and the group should be composed of two
        for match in self.cl[3:]:
            mol = readstring("smi",match)
            mol.addh()
            results = smarts.findall(mol)
            self.assertEquals(len(results), 2) # two carbon atoms have an L in both of these.
            for result in results:
                self.assertEquals(len(result), 2) #and each should have two atoms
        for not_match in self.not_cl:
            mol = readstring("smi",not_match)
            mol.addh()
            self.assertEquals(len(smarts.findall(mol)), 0)

if __name__ == "__main__":
    unittest.main()
