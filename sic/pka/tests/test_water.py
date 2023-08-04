"""
Tests whether the different protonation states of water have their pKa correctly identified.
"""

import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.structure import connectivity_table

class WaterPkaTest(unittest.TestCase):
    def setUp(self):
        self.waters = ["[OH-]","O","[OH3+]"]
    
    def testWaterPka(self):
        water_mols = []
        for w_string in self.waters:
            mol = readstring("smi",w_string)
            mol.addh()
            mol.connectivity_table = connectivity_table.ConnectivityTable(mol)
            pka.get_all_pka(mol)
            water_mols.append(mol)
        #get_all_pka modifies state so now let's iterate through atoms
        for atom in water_mols[0]:
            if atom.atomicnum == 1:
                self.assertTrue(pka.get_pka(atom.idx,water_mols[0]) == 52)
            else:
                self.assertTrue(pka.get_pka(atom.idx,water_mols[0]) == 15.7)
        for atom in water_mols[1]:
            if atom.atomicnum == 1:
                self.assertTrue(pka.get_pka(atom.idx,water_mols[1]) == 15.7)
            else:
                self.assertTrue(pka.get_pka(atom.idx,water_mols[1]) == -1.7)
        for atom in water_mols[2]:
            if atom.atomicnum == 1:
                self.assertTrue(pka.get_pka(atom.idx,water_mols[2])== -1.7)
            else:
                self.assertTrue(pka.get_pka(atom.idx,water_mols[2])== None)
def main():
    unittest.main()

if __name__ == "__main__":
    main()
