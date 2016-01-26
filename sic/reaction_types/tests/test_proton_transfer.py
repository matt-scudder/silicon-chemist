#!/usr/bin/python
#coding=utf-8
"""
Tests all the characteristics of a proton transfer:
    1. Whether the cross-check score is generated accurately for each of the four "tiers" of it (see proton_transfer.py)
    2. Whether rearrangement results in the correct structure
"""

from .. import proton_transfer
import segmentation.segmentation as segmentation
import structure.similarity as similarity
import pka.pka as pka
import unittest
from pybel import readstring,Smarts

class ProtonTransferTest(unittest.TestCase):
    def setUp(self):
        #I- trying to remove H from HF - ΔpKa = -13.2
        self.really_bad_reaction = readstring("smi","[I-].F")
        self.really_bad_reaction.addh()
        pka.get_all_pka(self.really_bad_reaction)
        self.really_bad_reaction_products = readstring("smi","I.[F-]")
        self.really_bad_reaction_products.addh()
        self.really_bad_sources = segmentation.label_sources(self.really_bad_reaction)
        self.really_bad_sinks = segmentation.label_sinks(self.really_bad_reaction)
        #I- trying to remove H from HCl - ΔpKa = -3
        self.uphill_reaction = readstring("smi","[I-].Cl")
        self.uphill_reaction.addh()
        pka.get_all_pka(self.uphill_reaction)
        self.uphill_reaction_products = readstring("smi","I.[Cl-]")
        self.uphill_reaction_products.addh()
        self.uphill_sources = segmentation.label_sources(self.uphill_reaction)
        self.uphill_sinks = segmentation.label_sinks(self.uphill_reaction)
        #Cl- trying to remove H from HI - ΔpKa = +3
        self.downhill_reaction = readstring("smi","[Cl-].I")
        self.downhill_reaction.addh()
        pka.get_all_pka(self.downhill_reaction)
        self.downhill_reaction_products = readstring("smi","Cl.[I-]")
        self.downhill_reaction_products.addh()
        self.downhill_sources = segmentation.label_sources(self.downhill_reaction)
        self.downhill_sinks = segmentation.label_sinks(self.downhill_reaction)
        #water taking H off HI - ΔpKa = 25.7 
        self.very_downhill_reaction = readstring("smi","[OH-].I")
        self.very_downhill_reaction.addh()
        pka.get_all_pka(self.very_downhill_reaction)
        self.very_downhill_reaction_products = readstring("smi","O.[I-]")
        self.very_downhill_reaction_products.addh()
        self.very_downhill_sources = segmentation.label_sources(self.very_downhill_reaction)
        self.very_downhill_sinks = segmentation.label_sinks(self.very_downhill_reaction)
        

    def testReallyBad(self):
        reaction = proton_transfer.ProtonTransfer([self.really_bad_sources[0]],[self.really_bad_sinks[0]])
        self.assertTrue(reaction.cross_check() == 0.0)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.really_bad_reaction,self.really_bad_reaction_products))

    def testUphill(self):
        reaction = proton_transfer.ProtonTransfer([self.uphill_sources[0]],[self.uphill_sinks[0]])
        self.assertTrue(reaction.cross_check() == 0.42)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.uphill_reaction,self.uphill_reaction_products))

    def testDownhill(self):
        reaction = proton_transfer.ProtonTransfer([self.downhill_sources[0]],[self.downhill_sinks[0]])
        self.assertTrue(reaction.cross_check() == 0.72)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.downhill_reaction,self.downhill_reaction_products))

    def testVeryDownhill(self):
        reaction = proton_transfer.ProtonTransfer([self.very_downhill_sources[0]],[self.very_downhill_sinks[0]])
        self.assertTrue(reaction.cross_check() == 1.0)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.very_downhill_reaction,self.very_downhill_reaction_products))

if __name__ == "__main__":
    unittest.main()
