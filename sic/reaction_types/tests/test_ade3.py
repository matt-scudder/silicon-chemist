import unittest

from openbabel.pybel import readstring

from sic.reaction_types import reaction_factory
from sic.segmentation import segmentation
from sic.structure import similarity, connectivity_table

class ADE3test(unittest.TestCase):
    """
    Checks weather "ADE3" reaction works on "C=C" by checking if we can convert "CC=C.I" to "CC(I)C" assuming that the "reaction.cross_check() = 0.3"
    Checks weather "ADE3" reaction works on "Z=C" by checking if we can convert "CC(C)=[O].Cl" to "CC(C)(Cl)O" assuming that the "reaction.cross_check() = 0.3"
    modify the reaction file to pick the right source for this test to work correctly "source = self.sources[?]""
    """

    def setUp(self):
        # C=C Addition
        self.double = readstring("smi","CC=C.I")
        self.double.addh()
        self.double.connectivity_table = connectivity_table.ConnectivityTable(self.double)
        self.double_sources = segmentation.label_sources(self.double)
        self.double_sinks = segmentation.label_sinks(self.double)
        self.double_products = readstring("smi","CC(I)C")
        self.double_products.addh()
        # Z=C Addition
        self.Z_double = readstring("smi","CC(C)=[O].Cl")
        self.Z_double.addh()
        self.Z_double.connectivity_table = connectivity_table.ConnectivityTable(self.Z_double)
        self.Z_double_sources = segmentation.label_sources(self.Z_double)
        self.Z_double_sinks = segmentation.label_sinks(self.Z_double)
        self.Z_double_products = readstring("smi","CC(C)(Cl)O")
        self.Z_double_products.addh()
             
    def testdouble(self):
        reaction = reaction_factory.produce_reaction("ADE3",self.double_sources,self.double_sinks)
        self.assertEqual(reaction.cross_check(), 0.3)
        reaction.rearrange()
        self.assertEqual(*similarity.normalize_mols([self.double,self.double_products]))
    
    def test_Z_double(self):
        reaction = reaction_factory.produce_reaction("ADE3",self.Z_double_sources,self.Z_double_sinks)
        self.assertEqual(reaction.cross_check(), 0.3)
        reaction.rearrange()
        error = "Inital Product = {}, Actual Product = {}".format(self.z_doubleBond_C, self.z_doubleBond_C_products)
        self.assertEqual(*similarity.normalize_mols([self.Z_double,self.Z_double_products]), error)
    
if __name__ == "__main__":
    unittest.main()

