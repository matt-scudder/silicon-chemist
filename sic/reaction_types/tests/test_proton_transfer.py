"""
Tests all the characteristics of a proton transfer:
    1. Whether the cross-check score is generated accurately for each of the four "tiers" of it (see proton_transfer.py)
    2. Whether rearrangement results in the correct structure
    3. Whether the proton transfer step works on "Z=C", Z=[O,S,N], as the first step of the "AE" reaction of "Z=C"
"""

# For the "test_ZdoubleBond" test, modify the reaction file to pick the right source for this test to work correctly "source = self.sources[?]""

import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.reaction_types import proton_transfer
from sic.segmentation import segmentation
from sic.structure import similarity
from sic.structure.connectivity_table import ConnectivityTable

class ProtonTransferTest(unittest.TestCase):
    def setUp(self):
        #I- trying to remove H from HF - ΔpKa = -13.2
        self.really_bad_reaction = readstring("smi","[I-].F")
        self.really_bad_reaction.addh()
        self.really_bad_reaction.connectivity_table = ConnectivityTable(self.really_bad_reaction)
        pka.get_all_pka(self.really_bad_reaction)
        self.really_bad_reaction_products = readstring("smi","I.[F-]")
        self.really_bad_reaction_products.addh()
        self.really_bad_sources = segmentation.label_sources(self.really_bad_reaction)
        self.really_bad_sinks = segmentation.label_sinks(self.really_bad_reaction)
        #I- trying to remove H from HCl - ΔpKa = -3
        self.uphill_reaction = readstring("smi","[I-].Cl")
        self.uphill_reaction.addh()
        self.uphill_reaction.connectivity_table = ConnectivityTable(self.uphill_reaction)
        pka.get_all_pka(self.uphill_reaction)
        self.uphill_reaction_products = readstring("smi","I.[Cl-]")
        self.uphill_reaction_products.addh()
        self.uphill_sources = segmentation.label_sources(self.uphill_reaction)
        self.uphill_sinks = segmentation.label_sinks(self.uphill_reaction)
        #Cl- trying to remove H from HI - ΔpKa = +3
        self.downhill_reaction = readstring("smi","[Cl-].I")
        self.downhill_reaction.addh()
        self.downhill_reaction.connectivity_table = ConnectivityTable(self.downhill_reaction)
        pka.get_all_pka(self.downhill_reaction)
        self.downhill_reaction_products = readstring("smi","Cl.[I-]")
        self.downhill_reaction_products.addh()
        self.downhill_sources = segmentation.label_sources(self.downhill_reaction)
        self.downhill_sinks = segmentation.label_sinks(self.downhill_reaction)
        #water taking H off HI - ΔpKa = 25.7 
        self.very_downhill_reaction = readstring("smi","[OH-].I")
        self.very_downhill_reaction.addh()
        self.very_downhill_reaction.connectivity_table = ConnectivityTable(self.very_downhill_reaction)
        pka.get_all_pka(self.very_downhill_reaction)
        self.very_downhill_reaction_products = readstring("smi","O.[I-]")
        self.very_downhill_reaction_products.addh()
        self.very_downhill_sources = segmentation.label_sources(self.very_downhill_reaction)
        self.very_downhill_sinks = segmentation.label_sinks(self.very_downhill_reaction)
        #CC(C)=O - trying to remove H from 
        self.z_doubleBond_C = readstring("smi","CC(C)=O.CC[OH2+]")
        self.z_doubleBond_C.addh()
        self.z_doubleBond_C.connectivity_table = ConnectivityTable(self.z_doubleBond_C)
        pka.get_all_pka(self.z_doubleBond_C)
        self.z_doubleBond_C_products = readstring("smi","C[C+](C)O.CCO")
        self.z_doubleBond_C_products.addh()
        self.z_doubleBond_C_sources = segmentation.label_sources(self.z_doubleBond_C)
        self.z_doubleBond_C_sinks = segmentation.label_sinks(self.z_doubleBond_C)
    
    def testReallyBad(self):
        reaction = proton_transfer.ProtonTransfer([self.really_bad_sources[0]],[self.really_bad_sinks[0]])
        self.assertEquals(reaction.cross_check(), 0.0)
        reaction.rearrange()
        self.assertEquals(*similarity.normalize_mols([self.really_bad_reaction,self.really_bad_reaction_products]))

    def testUphill(self):
        reaction = proton_transfer.ProtonTransfer([self.uphill_sources[0]],[self.uphill_sinks[0]])
        self.assertTrue(reaction.cross_check() > 0.20 and reaction.cross_check() < 0.22) #not gonna rely on equalities on math.exp()-derived floats...
        reaction.rearrange()
        self.assertEquals(*similarity.normalize_mols([self.uphill_reaction,self.uphill_reaction_products]))

    def testDownhill(self):
        reaction = proton_transfer.ProtonTransfer([self.downhill_sources[0]],[self.downhill_sinks[0]])
        self.assertTrue(reaction.cross_check() > 0.9 and reaction.cross_check() < 0.91)
        reaction.rearrange()
        self.assertEquals(*similarity.normalize_mols([self.downhill_reaction,self.downhill_reaction_products]))

    def testVeryDownhill(self):
        reaction = proton_transfer.ProtonTransfer([self.very_downhill_sources[0]],[self.very_downhill_sinks[0]])
        self.assertEquals(reaction.cross_check(), 1.0)
        reaction.rearrange()
        self.assertEquals(*similarity.normalize_mols([self.very_downhill_reaction,self.very_downhill_reaction_products]))
       
    def test_ZdoubleBond(self):
        reaction = proton_transfer.ProtonTransfer([self.z_doubleBond_C_sources[1]],[self.z_doubleBond_C_sinks[0]])
        self.assertGreater(reaction.cross_check(), -10)
        reaction.rearrange()
        print("product =",self.z_doubleBond_C)
        print("Actual product =", self.z_doubleBond_C_products)
        self.assertEquals(*similarity.normalize_mols([self.z_doubleBond_C,self.z_doubleBond_C_products]))

if __name__ == "__main__":
    unittest.main()
