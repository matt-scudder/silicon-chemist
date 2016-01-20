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
	log.debug(self.reactants.write("smiles"))
	Y = self.sources[0]["atoms"]["Y:"]
	H = self.sinks[0]["atoms"]["H"]
	L = self.sinks[0]["atoms"]["L"]
	print L["atom"].idx
	print Y["atom"].idx
	log.debug("Y: index: %s; H index %s" % (Y["atom"].idx,H["atom"].idx))
        struct_ops.make_bond(Y,H)
	log.debug(self.reactants.write("smiles"))
	log.debug("L index: %s; H index %s" % (L["atom"].idx,H["atom"].idx))
        struct_ops.break_bond(L,H)
	log.debug(self.reactants.write("smiles"))
	self.assertTrue(self.reactants.write("smiles") == self.products.write("smiles"))

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("arf").setLevel(logging.DEBUG)
    unittest.main()
