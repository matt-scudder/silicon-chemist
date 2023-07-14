import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.reaction_types import reaction_factory
from sic.segmentation import segmentation
from sic.structure import similarity, connectivity_table

class ANTest(unittest.TestCase):
    """
    Tests whether AN reactions work by checking whether a carbocation reaction adds in.
    This test will be updated further when rearrangements ar ein.
    """

    def setUp(self):
        self.cation = readstring("smi","C[C+](C)C.[Cl-]")
        self.cation.addh()
        pka.get_all_pka(self.cation)
        self.cation.connectivity_table = connectivity_table.ConnectivityTable(self.cation)
        self.cation_sources = segmentation.label_sources(self.cation)
        self.cation_sinks = segmentation.label_sinks(self.cation)
        self.cation_products = readstring("smi","CC(C)(C)Cl")
        self.cation_products.addh()

    def testCation(self):
        reaction = reaction_factory.produce_reaction("AN",self.cation_sources,self.cation_sinks)
        self.assertTrue(reaction.cross_check() == 1.0)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.cation,self.cation_products))

if __name__ == "__main__":
    unittest.main()
