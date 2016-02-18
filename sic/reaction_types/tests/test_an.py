#!/usr/bin/python
#coding=utf-8
import unittest
import pka.pka as pka
import reaction_types.reaction_type_factory as reaction_type_factory
import structure.similarity as similarity
import structure.struct_ops as struct_ops
import segmentation.segmentation as segmentation
from pybel import readstring

class ANTest(unittest.TestCase):
    """
    Tests whether AN reactions work by checking whether a carbocation reaction adds in.
    This test will be updated further when rearrangements ar ein.
    """

    def setUp(self):
        self.cation = readstring("smi","C[C+](C)C.[Cl-]")
        self.cation.addh()
        pka.get_all_pka(self.cation)
        struct_ops.generate_connectivity_table(self.cation)
        self.cation_sources = segmentation.label_sources(self.cation)
        self.cation_sinks = segmentation.label_sinks(self.cation)
        self.cation_products = readstring("smi","CC(C)(C)Cl")
        self.cation_products.addh()

    def testCation(self):
        reaction = reaction_type_factory.produce_reaction_type("AN",self.cation_sources,self.cation_sinks)
        self.assertTrue(reaction.cross_check() == 1.0)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.cation,self.cation_products))

if __name__ == "__main__":
    unittest.main()
