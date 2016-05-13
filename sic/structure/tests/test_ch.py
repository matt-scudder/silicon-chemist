#!/usr/bin/python
#coding=utf-8
"""
Tests whether the adjacent CH checks for eliminations actually work, as well as the get_H_bonds function.
"""
import structure.properties as properties
from pybel import readstring
import structure.connectivity_table as connectivity_table
import segmentation.segmentation as segmentation
import unittest

class CHTest(unittest.TestCase):
    def setUp(self):
        ch_check_mol = readstring("smi","CC(C)(C)Br.C[O-]")
        ch_check_mol.addh()
        ctable = connectivity_table.ConnectivityTable(ch_check_mol)
        ch_check_mol.connectivity_table = ctable
        self.ch_check_mol = ch_check_mol
        self.sinks = segmentation.label_sinks(ch_check_mol)

    def testGetHBonded(self):
        self.assertTrue(len(properties.get_H_bonds(1,self.ch_check_mol)) == 3)
        self.assertTrue(len(properties.get_H_bonds(5,self.ch_check_mol)) < 1)

    def testCH(self):
        #first get the CL and Y
        actual_sink = False
        for sink in self.sinks:
            if sink.subtype == "C-L":
                actual_sink = sink
                break
        #ok, now check CH
        CH = properties.get_adjacent_ch(actual_sink)
        #make sure it is actually a CH
        self.assertTrue(CH) #since atom indices start at 1, CH should always be truthy.
        self.assertTrue(len(properties.get_H_bonds(CH,self.ch_check_mol)) == 3)
        #check against Br now
        possible_L_H = properties.get_adjacent_ch(actual_sink,carbon_label="L") #tests the flexibility of our function, too
        self.assertFalse(possible_L_H)

if __name__ == "__main__":
    unittest.main()
