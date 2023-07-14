import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.reaction_types import reaction_factory
from sic.segmentation import segmentation
from sic.structure import similarity, connectivity_table

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
        self.primary.connectivity_table = connectivity_table.ConnectivityTable(self.primary)
        self.primary_sources = segmentation.label_sources(self.primary)
        self.primary_sinks = segmentation.label_sinks(self.primary)
        self.primary_products = readstring("smi","CCO.[Cl-]")
        self.primary_products.addh()
        #secondary
        self.secondary = readstring("smi","CC(C)Cl.[OH-]")
        self.secondary.addh()
        pka.get_all_pka(self.secondary)
        self.secondary.connectivity_table = connectivity_table.ConnectivityTable(self.secondary)
        self.secondary_sources = segmentation.label_sources(self.secondary)
        self.secondary_sinks = segmentation.label_sinks(self.secondary)
        self.secondary_products = readstring("smi","CC(C)O.[Cl-]")
        self.secondary_products.addh()
        #tertiary
        self.tertiary = readstring("smi","CC(C)(C)Cl.[OH-]")
        self.tertiary.addh()
        pka.get_all_pka(self.tertiary)
        self.tertiary.connectivity_table = connectivity_table.ConnectivityTable(self.tertiary)
        self.tertiary_sources = segmentation.label_sources(self.tertiary)
        self.tertiary_sinks = segmentation.label_sinks(self.tertiary)
        #no products because this one shouldn't happen.

    #TODO: Update the below when we get C-H...
    def testPrimary(self):
        reaction = reaction_factory.produce_reaction("SN2",self.primary_sources,self.primary_sinks)
        self.assertTrue(reaction.cross_check() == 0.95)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.primary,self.primary_products))

    def testSecondary(self):
        reaction = reaction_factory.produce_reaction("SN2",self.secondary_sources,self.secondary_sinks)
        self.assertTrue(reaction.cross_check() == 0.6)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.secondary,self.secondary_products))

    def testTertiary(self):
        reaction = reaction_factory.produce_reaction("SN2",self.tertiary_sources,self.tertiary_sinks)
        self.assertTrue(reaction.cross_check() == 0.0)
        

if __name__ == "__main__":
    unittest.main()
