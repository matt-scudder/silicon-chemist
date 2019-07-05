import unittest
import pybel
from pybel import readstring,Smarts
import structure.similarity as similarity
import structure.properties as properties
import structure.connectivity_table as connectivity_table
import segmentation.segmentation as segmentation
import segmentation.sources as sources
import segmentation.source as source
import reaction_types.reaction_factory as reaction_factory
# Checks weather "ADE3" reaction works:
#  if we can convert "CC=C.IBr" to "C1C[I+]1.[Br-]" to test the ("C=C", "Y-L") interaction
# Modify the reaction file to pick the right source and sink for this test to work correctly "source = self.sources[?] "
class NuLTest(unittest.TestCase):
    def setUp(self):
        self.double_Y_L= readstring("smi","C=C.IBr")
        self.double_Y_L.addh()
        self.double_Y_L.connectivity_table = connectivity_table.ConnectivityTable(self.double_Y_L)
        self.double_Y_L_sources = segmentation.label_sources(self.double_Y_L)
        print "Sources = ",self.double_Y_L_sources
        self.double_Y_L_sinks = segmentation.label_sinks(self.double_Y_L)
        print "Sinks =", self.double_Y_L_sinks
        self.double_Y_L_products = readstring("smi","C1C[I+]1.[Br-]")
        self.double_Y_L_products.addh()

    def testdouble(self):
        reaction = reaction_factory.produce_reaction("NuL",self.double_Y_L_sources,self.double_Y_L_sinks)
        print "reaction.cross_check() =", reaction.cross_check()
        self.assertTrue(reaction.cross_check() == 0.3)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.double_Y_L,self.double_Y_L_products))

if __name__ == "__main__":
    unittest.main()