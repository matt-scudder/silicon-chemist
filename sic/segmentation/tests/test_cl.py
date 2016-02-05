#!/usr/bin/python
#coding=utf-8
"""
Tests whether C-L can be identified.
"""

import segmentation.sinks as sinks
from pybel import readstring,Smarts
import unittest

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
            self.assertTrue(len(results) == 1) #only one set of atoms 
        for match in self.cl[3:]:
            mol = readstring("smi",match)
            mol.addh()
            results = smarts.findall(mol)
            self.assertTrue(len(results) == 2) # two carbon atoms have an L in both of these.
        for not_match in self.not_cl:
            mol = readstring("smi",not_match)
            mol.addh()
            self.assertTrue(len(smarts.findall(mol)) == 0)

if __name__ == "__main__":
    unittest.main()
