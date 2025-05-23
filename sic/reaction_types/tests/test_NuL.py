import unittest

from openbabel.pybel import readstring

from sic.reaction_types import reaction_factory
from sic.segmentation import segmentation
from sic.structure import similarity, connectivity_table

class NuLTest(unittest.TestCase):
    """
    Checks weather "ADE3" reaction works:
      if we can convert "CC=C.IBr" to "C1C[I+]1.[Br-]" to test the ("C=C", "Y-L") interaction
    Modify the reaction file to pick the right source and sink for this test to work correctly "source = self.sources[?] "
    """

    def setUp(self):
        self.double_Y_L= readstring("smi","C=C.IBr")
        self.double_Y_L.addh()
        self.double_Y_L.connectivity_table = connectivity_table.ConnectivityTable(self.double_Y_L)
        self.double_Y_L_sources = segmentation.label_sources(self.double_Y_L)
        print("Sources = ",self.double_Y_L_sources)
        self.double_Y_L_sinks = segmentation.label_sinks(self.double_Y_L)
        print("Sinks =", self.double_Y_L_sinks)
        self.double_Y_L_products = readstring("smi","C1C[I+]1.[Br-]")
        self.double_Y_L_products.addh()

    def testdouble(self):
        reaction = reaction_factory.produce_reaction("NuL",self.double_Y_L_sources,self.double_Y_L_sinks)
        print("reaction.cross_check() =", reaction.cross_check())
        self.assertEqual(reaction.cross_check(), 0.3)
        reaction.rearrange()
        self.assertEqual(*similarity.normalize_mols([self.double_Y_L,self.double_Y_L_products]))

if __name__ == "__main__":
    unittest.main()