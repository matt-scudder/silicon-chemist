#!/usr/bin/python
#coding=utf-8
"""
Tests whether the make_bond and break_bond methods work for reactions that move a single atom,
e.g. proton transfers.
Because these methods are crucial to our modification of the reaction state,
it is very important that they work correctly.

This module also tests the label_sources and label_sinks methods in segmentation,
as they are necessary to get everything into the right format, and the tanimoto coefficient.
"""

from .. import struct_ops, similarity
from ...segmentation import segmentation 
import unittest
from pybel import readstring,Smarts
import logging
import sys

class SingleAtomMovementTest(unittest.TestCase):
    def setUp(self):
        acid = readstring("smi","F")
        acid.addh()
        self.acid = acid
        base = readstring("smi","O")
        base.addh()
        self.base = base
        conj_base = readstring("smi","[F-]")
        conj_base.addh()
        self.conj_base = conj_base
        conj_acid = readstring("smi","[OH3+]")
        self.conj_acid = conj_acid
        #since we know that acid is sink and base is source, let's just do normal labeling
        self.sources = segmentation.label_sources(base)
        self.sinks = segmentation.label_sinks(acid)

    def testSingleAtomMovement(self):
        #this test has two parts, because we *must* ensure that bond breaking happens before bond breaking
        #but tests are normally run in any order.
        log = logging.getLogger("arf")
        struct_ops.make_bond(self.sources[0]["atoms"]["Y:"],self.sinks[0]["atoms"]["H"],single_atom=True)
        log.debug(self.base.write("smiles"))
        log.debug(self.conj_acid.write("smiles"))
        self.assertTrue(similarity.tanimoto(self.base,self.conj_acid) == 1.0) #because self.base got modified...
        struct_ops.break_bond(self.sinks[0]["atoms"]["L"],self.sinks[0]["atoms"]["H"],single_atom=True)
        log.debug(self.acid.write("smiles"))
        self.assertTrue(similarity.tanimoto(self.acid,self.conj_base) == 1.0)

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("arf").setLevel(logging.DEBUG)
    unittest.main()
