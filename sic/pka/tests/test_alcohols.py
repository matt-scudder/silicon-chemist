"""
Tests whether alcohols are correctly identified, and that the degrees are not confused with each other.
"""

import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.structure import connectivity_table

class AlcoholPkaTest(unittest.TestCase):
    def setUp(self):
        self.alcohols = {
            "CO": "methyl",
            "CCO": "primary",
            "CC(C)O": "secondary",
            "CC(O)(C)C": "tert",
        }
    
    def testAlcoholPka(self):
        for smiles_string, alcohol_name in self.alcohols.items():
            with self.subTest(alcohol_name):
                mol = readstring("smi",smiles_string)
                mol.addh()
                mol.connectivity_table = connectivity_table.ConnectivityTable(mol)
                pka.get_all_pka(mol)

                #get_all_pka modifies state so now let's iterate through atoms
                match alcohol_name:
                    case "methyl":
                        for atom in mol:
                            if atom.idx == 6: #because doing the connectivity test is a pain...
                                self.assertEquals(pka.get_pka(atom.idx,mol), 15.5)
                            elif atom.atomicnum == 8:
                                self.assertEquals(pka.get_pka(atom.idx,mol), -2.4)
                    case "primary":
                        for atom in mol:
                            if atom.idx == 9:
                                self.assertEquals(pka.get_pka(atom.idx,mol), 16)
                            elif atom.atomicnum == 8:
                                self.assertEquals(pka.get_pka(atom.idx,mol), -2.4)
                    case "secondary":
                        for atom in mol:
                            if atom.idx == 12:
                                self.assertEquals(pka.get_pka(atom.idx,mol), 18)
                            elif atom.atomicnum == 8:
                                self.assertEquals(pka.get_pka(atom.idx,mol), -2.4)
                    case "tert":
                        for atom in mol:
                            if atom.idx == 9:
                                self.assertEquals(pka.get_pka(atom.idx,mol), 19)
                            elif atom.atomicnum == 8:
                                self.assertEquals(pka.get_pka(atom.idx,mol), -2.4)
def main():
    unittest.main()

if __name__ == "__main__":
    main()
