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

import structure.struct_ops as struct_ops
import structure.similarity as similarity
import segmentation.segmentation as segmentation
import unittest
from pybel import readstring,Smarts
import logging
import sys

#NOTE: impossibility of ΔpKa < -10 reactions is left to the logic that calls this code...
class SingleAtomMovementTest(unittest.TestCase):
    def setUp(self):
        reactants = readstring("smi","F.O")
        reactants.addh()
        self.reactants = reactants
        products = readstring("smi","[F-].[OH3+]")
        products.addh()
        self.products = products
        self.sources = segmentation.label_sources(reactants)
        self.sinks = segmentation.label_sinks(reactants)

    def testSingleAtomMovement(self):
        #this test has two parts, because we *must* ensure that bond breaking happens before bond breaking
        #but tests are normally run in any order.
        log = logging.getLogger("arf")
        log.debug(self.reactants.write("can"))
        Y = self.sources[0].get_atom("Y")
        H = self.sinks[0].get_atom("H")
        L = self.sinks[0].get_atom("L")
        log.debug("Y index: %s; H index %s" % (Y,H))
        struct_ops.make_bond(Y,H,self.reactants)
        log.debug(self.reactants.write("can"))
        log.debug("L index: %s; H index %s" % (L,H))
        struct_ops.break_bond(H,L,self.reactants)
        log.debug(self.reactants.write("can"))
        self.assertTrue(similarity.is_same_molecule(self.reactants,self.products))

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("arf").setLevel(logging.DEBUG)
    unittest.main()
