"""
This unit tests whether non-halogen lone pairs are identified correctly by our source patterns.
Because we're testing the PATTERNS here rather than the segmentation functions themselves,
we don't call the functions in segmentation.py.
"""

import unittest

from openbabel.pybel import readstring,Smarts

from sic.segmentation import sources

class LonePairIdentificationTest(unittest.TestCase):
    def setUp(self):
        self.SOURCES = sources.SOURCES
        self.oxygen_pairs = ["O","[OH-]","CO","CC=O","COC"] #other connectivities are similar enough to these to not need testing
        self.full_oxygen_valence = ["[OH3+]","C[OH2+]","CC=[OH+]","C[OH+]C"]
        self.nitrogen_pairs = ["N","CN","CNC","CN(C)C"] #all have nitrogen lone pair, other tests unnecessary
        self.full_nitrogen_valence = ["[NH4+]","C[NH3+]","C[NH2+]C","C[NH+](C)C","C[N+](C)(C)C"]
        self.sulfur_pairs = ["CCS","CSC"] #only ever see thioesters and thiols in organic
        self.full_sulfur_valence = ["CC[SH2+]","C[SH+]C"]

    def testOxygenSources(self):
        """
        Tests whether we can identify oxygen sources without also identifying oxygens with no lone pairs.
        """

        smarts = Smarts(self.SOURCES["Y"])
        for oxygen_pair in self.oxygen_pairs:
            mol = readstring("smi",oxygen_pair)
            mol.addh()
            results = smarts.findall(mol)
            self.assertEquals(len(results), 1) #only one atom that's relevant
            self.assertEquals(len(results[0]), 1) #only catch the lone pair atom

        for full_valence in self.full_oxygen_valence:
            mol = readstring("smi",full_valence)
            mol.addh()
            self.assertEquals(len(smarts.findall(mol)), 0) #full valence means you shouldn't get anything
    
    def testNitrogenSources(self):
        """
        Tests whether we can identify nitrogen sources without also identifying nitrogens with no lone pairs.
        """

        smarts = Smarts(self.SOURCES["Y"])
        for nitrogen_pair in self.nitrogen_pairs:
            mol = readstring("smi",nitrogen_pair)
            mol.addh()
            results = smarts.findall(mol)
            self.assertEquals(len(results), 1) #only one atom that's relevant
            self.assertEquals(len(results[0]), 1) #only catch the lone pair atom

        for full_valence in self.full_nitrogen_valence:
            mol = readstring("smi",full_valence)
            mol.addh()
            self.assertEquals(len(smarts.findall(mol)), 0) #full valence means you shouldn't get anything

    def testSulfurSources(self):
        """
        Tests whether we can identify sulfur sources without also identifying sulfurs with no lone pairs.
        """

        smarts = Smarts(self.SOURCES["Y"])
        for sulfur_pair in self.sulfur_pairs:
            mol = readstring("smi",sulfur_pair)
            mol.addh()
            results = smarts.findall(mol)
            self.assertEquals(len(results), 1) #only one atom that's relevant
            self.assertEquals(len(results[0]), 1) #only catch the lone pair atom

        for full_valence in self.full_sulfur_valence:
            mol = readstring("smi",full_valence)
            mol.addh()
            self.assertEquals(len(smarts.findall(mol)), 0) #full valence means you shouldn't get anything

def main():
    unittest.main()

if __name__ == "__main__":
    main()