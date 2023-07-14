import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.reaction_types import reaction_factory
from sic.segmentation import segmentation
from sic.structure import similarity, connectivity_table

class E2Test(unittest.TestCase):
    """
    Tests whether E2 reactions work by checking whether we can eliminate from t-butyl chloride.
    """

    def setUp(self):
        self.t_butyl = readstring("smi","CC(C)(C)Cl.C[O-]")
        self.t_butyl.addh()
        pka.get_all_pka(self.t_butyl)
        self.t_butyl.connectivity_table = connectivity_table.ConnectivityTable(self.t_butyl)
        self.t_butyl_sources = segmentation.label_sources(self.t_butyl)
        self.t_butyl_sinks = segmentation.label_sinks(self.t_butyl)
        self.t_butyl_products = readstring("smi","C=C(C)C.[Cl-].CO")
        self.t_butyl_products.addh()

    def testCation(self):
        reaction = reaction_factory.produce_reaction("E2",self.t_butyl_sources,self.t_butyl_sinks)
        self.assertTrue(reaction.cross_check() == 0.3)
        reaction.rearrange()
        #print self.t_butyl.write("can")
        self.assertTrue(similarity.is_same_molecule(self.t_butyl,self.t_butyl_products))

if __name__ == "__main__":
    unittest.main()
