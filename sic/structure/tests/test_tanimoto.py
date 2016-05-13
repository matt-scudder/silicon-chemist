#!/usr/bin/python
#coding=utf-8
"""
Tests whether the Tanimoto coefficient is an effective way of figuring out whether we have the product.
Tests for both the individual case, and for the "reaction state" case.

Since the properties of a ReactionState object that aren't the list of reactants are not relevant
to this test, we will only use a list of reactants.
"""

from .. import similarity
import unittest
from pybel import readstring

class TanimotoSimilarityTest(unittest.TestCase):
    def setUp(self):
        self.complicated_molecule = readstring("smi","OCCC=CCCCCCC(O)CCCCCCCCC(O)CO") #long enough that 7-atom linear may fail
        self.complicated_molecule.addh() #just in case - remember we always addh so we need to make sure it works with addh
        self.slightly_different_complicated_molecule = readstring("smi","OCCC=CCCCCCC([O-])CCCCCCCCC(O)CO")
        self.slightly_different_complicated_molecule.addh()
        self.identical_complicated_molecule = readstring("smi","OCCC=CCCCCCC(O)CCCCCCCCC(O)CO")
        self.identical_complicated_molecule.addh()
        base_react_state = []
        same_react_state = []
        different_react_state = []
        mol1 = readstring("smi","OCC(CO)C[O-]")
        mol1.addh()
        base_react_state.append(mol1)
        mol2 = readstring("smi","[F-]")
        mol2.addh()
        base_react_state.append(mol2)
        mol3 = readstring("smi","CC(=O)[O-]")
        mol3.addh()
        base_react_state.append(mol3)
        same_react_state.append(mol3)
        same_react_state.append(mol1)
        same_react_state.append(mol2)
        different_react_state.append(mol1)
        different_react_state.append(mol2)
        mol4 = readstring("smi","CC(=O)O")
        mol4.addh()
        different_react_state.append(mol4)
        self.base_react_state = base_react_state
        self.same_react_state = same_react_state
        self.different_react_state = different_react_state
        self.mol1 = mol1
        self.mol2 = mol2
        self.mol3 = mol3
        self.mol4 = mol4


    def testComplicatedMolecule(self):
        self.assertTrue(similarity.tanimoto(self.complicated_molecule,self.slightly_different_complicated_molecule) != 1.0)
        self.assertTrue(similarity.tanimoto(self.complicated_molecule,self.identical_complicated_molecule) == 1.0)
        print similarity.tanimoto(self.mol3,self.mol3)

    def testTanimotoMapping(self):
        self.assertTrue(similarity.is_same_state(self.base_react_state,self.same_react_state))
        self.assertFalse(similarity.is_same_state(self.base_react_state,self.different_react_state))

def main():
    unittest.main()

if __name__ == "__main__":
    main()

