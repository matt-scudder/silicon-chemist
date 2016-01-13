#!/usr/bin/python
#coding=utf-8
"""
This unit test tests whether the halogen acids and their conjugate bases
are identified properly as sources and sinks.
Because we're testing the PATTERNS here rather than the segmentation functions themselves,
we don't call the functions in segmentation.py.
"""
from .. import sources, sinks
from pybel import readstring,Smarts
import unittest

class HalogenIdentificationTest(unittest.TestCase):
    def setUp(self):
        self.halogen_c_bases = ["[Cl-]","[F-]","[Br-]","[I-]"]
        self.halogen_acids = ["Cl","F","Br","I"]
        self.halogen_not_match = ["CCl","CF","CBr","CI"] #these shouldn't match either Y: or H-L
        self.SOURCES = sources.SOURCES
        self.SINKS = sinks.SINKS
        
    def testHalogenSources(self):
        """
        Tests whether the haloges are correctly identified as sources when they are,
        and correctly NOT identified as sources when they aren't.
        """
        smarts = Smarts(self.SOURCES["Y:"])
        for c_base in self.halogen_c_bases:
            mol = readstring("smi",c_base)
            mol.addh()
            results = smarts.findall(mol)
            self.assertTrue(len(results) == 1) #should be only one group that's a source
            self.assertTrue(len(results[0]) == 1) #and only 1 atom in the group
        
        for acid in self.halogen_acids:
            mol = readstring("smi",acid)
            mol.addh()
            self.assertTrue(len(smarts.findall(mol)) == 0) #shouldn't find anything here

        for not_match in self.halogen_not_match:
            mol = readstring("smi",not_match)
            mol.addh()
            self.assertTrue(len(smarts.findall(mol)) == 0) #shouldn't find anything here

    def testHalogenSinks(self):
        """
        Same as for testHalogenSources, but this time for H-L sinks.
        """
        smarts = Smarts(self.SINKS["H-L"])
        for c_base in self.halogen_c_bases:
            mol = readstring("smi",c_base)
            mol.addh()
            self.assertTrue(len(smarts.findall(mol)) == 0) #shouldn't find anything here

        for acid in self.halogen_acids:
            mol = readstring("smi",acid)
            mol.addh()
            results = smarts.findall(mol)
            self.assertTrue(len(results) == 1) #should be exactly 1 acid
            self.assertTrue(len(results[0]) == 2) #catch the H for sure

        for not_match in self.halogen_not_match:
            mol = readstring("smi",not_match)
            mol.addh()
            self.assertTrue(len(smarts.findall(mol)) == 0) #shouldn't find anything here

def main():
    unittest.main()

if __name__ == "__main__":
    main()
