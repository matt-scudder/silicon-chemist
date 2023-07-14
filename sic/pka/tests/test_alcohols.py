"""
Tests whether alcohols are correctly identified, and that the degrees are not confused with each other.
"""

import unittest
from openbabel.pybel import readstring
from .. import pka_chart,pka

class AlcoholPkaTest(unittest.TestCase):
    def setUp(self):
        self.PKA_CHART = pka_chart.PKA_CHART
        self.alcohols = ["CO","CCO","CC(C)O","CC(O)(C)C"] #methyl, primary, secondary, tert
    
    def testAlcoholPka(self):
        alcohol_mols = []
        for w_string in self.alcohols:
            mol = readstring("smi",w_string)
            mol.addh()
            pka.get_all_pka(mol)
            alcohol_mols.append(mol)
        #get_all_pka modifies state so now let's iterate through atoms
        for atom in alcohol_mols[0]:
            if atom.idx == 6: #because doing the connectivity test is a pain...
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[0]) == 15.5)
            elif atom.atomicnum == 8:
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[0]) == -2.4)
        for atom in alcohol_mols[1]:
            if atom.idx == 9:
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[1]) == 16)
            elif atom.atomicnum == 8:
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[1]) == -2.4)
        for atom in alcohol_mols[2]:
            if atom.idx == 12:
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[2])== 18)
            elif atom.atomicnum == 8:
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[2])== -2.4)
        for atom in alcohol_mols[3]:
            if atom.idx == 9:
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[3])== 19)
            elif atom.atomicnum == 8:
                self.assertTrue(pka.get_pka(atom.idx,alcohol_mols[3])== -2.4)
def main():
    unittest.main()

if __name__ == "__main__":
    main()
