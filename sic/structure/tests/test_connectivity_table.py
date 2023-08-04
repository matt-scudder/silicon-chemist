
import unittest

from openbabel import pybel

from sic.segmentation import segmentation
from sic.structure import struct_ops, connectivity_table

class ConnectivityTableTest(unittest.TestCase):
    """
    Tests the connectivity table functions for use with "closer to product" checking.
    """

    def setUp(self):
        mol = pybel.readstring("smi","O.F")
        mol.addh()
        ctable = connectivity_table.ConnectivityTable(mol)
        mol.connectivity_table = ctable
        self.sources = segmentation.label_sources(mol)
        self.sinks = segmentation.label_sinks(mol)
        self.mol = mol
    
    def testTable(self):
        """
        Tests whether the table is generated correctly at all.
        """

        ctable = self.mol.connectivity_table #typing is hard
        #F bond
        
        self.assertIn(5, ctable.get_atoms_bonded(2))
        self.assertIn(2, ctable.get_atoms_bonded(5))
        #O bond
        self.assertIn(3, ctable.get_atoms_bonded(1))
        self.assertIn(1, ctable.get_atoms_bonded(3))

    def testAddBond(self):
        """
        Tests whether adding a bond updates the table. Note that readstring is deterministic, so it
        doesn't matter whether we add or remove first, the indices are the same.
        """

        ctable = self.mol.connectivity_table #typing is hard
        #sinks[2] is the F, not the O.
        O = self.sources[0].get_atom("Y")
        H = self.sinks[2].get_atom("H")
        struct_ops.make_bond(O,H,self.mol)
        self.assertIn(5, ctable.get_atoms_bonded(1))
        self.assertIn(1, ctable.get_atoms_bonded(5))

    def testRemoveBond(self):
        """
        Tests whether removing a bond updates the table.
        """

        ctable = self.mol.connectivity_table #typing is hard
        F = self.sinks[2].get_atom("L")
        H = self.sinks[2].get_atom("H")
        struct_ops.break_bond(H,F,self.mol)
        self.assertNotIn(5, ctable.get_atoms_bonded(2))
        self.assertNotIn(2, ctable.get_atoms_bonded(5))

if __name__ == "__main__":
    unittest.main()
