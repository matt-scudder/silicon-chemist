#!/usr/bin/python
#coding=utf-8

import unittest
import structure.struct_ops as struct_ops
import structure.connectivity_table as connectivity_table
import pybel
import segmentation.segmentation as segmentation

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
        self.assertTrue(2 in ctable and 5 in ctable[2])
        self.assertTrue(5 in ctable and 2 in ctable[5])
        #O bond
        self.assertTrue(1 in ctable and 3 in ctable[1])
        self.assertTrue(3 in ctable and 1 in ctable[3])

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
        self.assertTrue(1 in ctable and 5 in ctable[1])
        self.assertTrue(5 in ctable and 1 in ctable[5])

    def testRemoveBond(self):
        """
        Tests whether removing a bond updates the table.
        """
        ctable = self.mol.connectivity_table #typing is hard
        F = self.sinks[2].get_atom("L")
        H = self.sinks[2].get_atom("H")
        struct_ops.break_bond(H,F,self.mol)
        self.assertFalse(2 in ctable and 5 in ctable[2])
        self.assertFalse(5 in ctable and 2 in ctable[5])

if __name__ == "__main__":
    unittest.main()
