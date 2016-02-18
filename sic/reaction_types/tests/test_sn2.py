#!/usr/bin/python
#coding=utf-8
import unittest
import pka.pka as pka
import reaction_types.reaction_type_factory as reaction_type_factory
import structure.similarity as similarity
import structure.struct_ops as struct_ops
import segmentation.segmentation as segmentation
from pybel import readstring

class SN2Test(unittest.TestCase):
    """
    Tests whether SN2 works by checking whether a primary and secondary SN2 function
    and provide correct scores, and whether a tertiary SN2 gives a 0.0.
    """

    def setUp(self):
        #three reactions - one primary, one secondary, one tertiary
        #same Nu, same L, different grease
        #primary, ΔpKa >> 10
        self.primary = readstring("smi","CCCl.[OH-]")
        self.primary.addh()
        pka.get_all_pka(self.primary)
        struct_ops.generate_connectivity_table(self.primary)
        self.primary_sources = segmentation.label_sources(self.primary)
        self.primary_sinks = segmentation.label_sinks(self.primary)
        self.primary_products = readstring("smi","CCO.[Cl-]")
        self.primary_products.addh()
        #secondary
        self.secondary = readstring("smi","CC(C)Cl.[OH-]")
        self.secondary.addh()
        pka.get_all_pka(self.secondary)
        struct_ops.generate_connectivity_table(self.secondary)
        self.secondary_sources = segmentation.label_sources(self.secondary)
        self.secondary_sinks = segmentation.label_sinks(self.secondary)
        self.secondary_products = readstring("smi","CC(C)O.[Cl-]")
        self.secondary_products.addh()
        #tertiary
        self.tertiary = readstring("smi","CC(C)(C)Cl.[OH-]")
        self.tertiary.addh()
        pka.get_all_pka(self.tertiary)
        struct_ops.generate_connectivity_table(self.tertiary)
        self.tertiary_sources = segmentation.label_sources(self.tertiary)
        self.tertiary_sinks = segmentation.label_sinks(self.tertiary)
        #no products because this one shouldn't happen.

    #TODO: Update the below when we get C-H...
    def testPrimary(self):
        reaction = reaction_type_factory.produce_reaction_type("SN2",self.primary_sources,self.primary_sinks)
        self.assertTrue(reaction.cross_check() == 1.0)
        reaction.rearrange()
        print self.primary.write("smiles")
        self.assertTrue(similarity.is_same_molecule(self.primary,self.primary_products))

    def testSecondary(self):
        reaction = reaction_type_factory.produce_reaction_type("SN2",self.secondary_sources,self.secondary_sinks)
        self.assertTrue(reaction.cross_check() == 0.6)
        reaction.rearrange()
        print self.secondary.write("smiles")
        self.assertTrue(similarity.is_same_molecule(self.secondary,self.secondary_products))

    def testTertiary(self):
        reaction = reaction_type_factory.produce_reaction_type("SN2",self.tertiary_sources,self.tertiary_sinks)
        self.assertTrue(reaction.cross_check() == 0.0)
        

if __name__ == "__main__":
    unittest.main()
