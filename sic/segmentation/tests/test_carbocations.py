"""
Tests whether carbocations can be identified.
"""

import unittest

from openbabel.pybel import readstring,Smarts

from sic.segmentation import sinks

class CarbocationIdentificationTest(unittest.TestCase):
    def setUp(self):
        self.carbocations = ["CCC[CH2+]","[CH3+]","C[CH+]C","C[C+](C)C"]
        self.not_carbocations = ["CCCC","C","CCC","CC(C)C"]
        self.SINKS = sinks.SINKS

    def testCarbocationMatch(self):
        """
        Tests whether carbocations are identified correctly as sinks when they are,
        and NOT when there's no carbocation.
        """

        smarts = Smarts(self.SINKS["C+"])
        for carbocation in self.carbocations:
            mol = readstring("smi",carbocation)
            mol.addh()
            results = smarts.findall(mol)
            self.assertEquals(len(results), 1) #only one molecule
            self.assertEquals(len(results[0]), 1) # only one carbon atom
        
        for not_match in self.not_carbocations:
            mol = readstring("smi",not_match)
            mol.addh()
            self.assertEquals(len(smarts.findall(mol)), 0)

if __name__ == "__main__":
    unittest.main()
