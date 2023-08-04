"""
Tests whether attaching two groups together works.
Uses a carbocation bonding to an O lone pair to test whether the two groups bond to each other.
O- is used so that there's only one source and one sink here.
"""

import unittest

from openbabel.pybel import readstring

from sic.structure import struct_ops, similarity, connectivity_table
from sic.segmentation import segmentation

class GroupAttachmentTest(unittest.TestCase):
    def setUp(self):
        reactants = readstring("smi","CCCCCC[O-].C[C+](C)C")
        products = readstring("smi","CCCCCCOC(C)(C)C")
        reactants.addh()
        products.addh()
        ctable = connectivity_table.ConnectivityTable(reactants)
        reactants.connectivity_table = ctable
        self.reactants = reactants
        self.products = products
        self.sources = segmentation.label_sources(reactants)
        self.sinks = segmentation.label_sinks(reactants)

    def testAttachment(self):
        #make bond from [O-] to the C+
        Y = self.sources[0].get_atom("Y")
        C = self.sinks[0].get_atom("C+")
        print(Y)
        print(C)
        print(self.reactants.write("can"))
        struct_ops.make_bond(Y,C,self.reactants) #self.reactants is a Molecule object
        print(self.reactants.write("can"))
        self.assertEquals(*similarity.normalize_mols([self.reactants,self.products]))

if __name__ == "__main__":
    unittest.main()

