import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.reaction_types import reaction_factory
from sic.segmentation import segmentation
from sic.structure import similarity, connectivity_table

class ADN(unittest.TestCase):
    """
    Checks weather "ADN" reaction works by checking if we can convert "N#C=C.CC[O-]" to "N#C[CH-]COCC"
    Differnt input format : "C=CC#N.CC[O-]"
    test "Z=C": CC=O.[CH2-]C=O  to CC([O-])CC=O
    Modify the reaction file to pick the right source and sink for this test to work correctly "source = self.sources[?] "
    """

    def setUp(self):
        # C=C Addition
        self.double = readstring("smi","C=OC#N.[O-]CC")
        self.double.addh()
        self.double.connectivity_table = connectivity_table.ConnectivityTable(self.double)
        pka.get_all_pka(self.double)
        self.double_sources = segmentation.label_sources(self.double)
        self.double_sinks = segmentation.label_sinks(self.double)
        self.double_products = readstring("smi","C(OCC)[O-]C#N")
        self.double_products.addh()
        #print "molecule =", utils.write_all_bonds(self.double)
        # Z=C Addition
        self.Z_double = readstring("smi","CC=O.[CH2-]C=O")
        self.Z_double.addh()
        self.Z_double.connectivity_table = connectivity_table.ConnectivityTable(self.Z_double)
        pka.get_all_pka(self.Z_double)
        self.Z_double_sources = segmentation.label_sources(self.Z_double)
        self.Z_double_sinks = segmentation.label_sinks(self.Z_double)
        self.Z_double_products = readstring("smi","CC([O-])CC=O")
        self.Z_double_products.addh()

    # tests "C=C" 
    def testdouble(self):
        reaction = reaction_factory.produce_reaction("ADN",self.double_sources,self.double_sinks)
        self.assertEquals(reaction.cross_check(), 1)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.double,self.double_products))
    
    #tests "Z=C"
    def test_Z_double(self):
        reaction = reaction_factory.produce_reaction("ADN",self.Z_double_sources,self.Z_double_sinks)
        self.assertLess(reaction.cross_check(), 1)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.Z_double,self.Z_double_products))

if __name__ == "__main__":
    unittest.main()
