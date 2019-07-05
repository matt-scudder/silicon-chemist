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
# Checks weather "ADE3" reaction works by checking these cases:
#1] if we can convert "CC=C.I" to "C[C+]C.[I-]" to test the ("C=C", "H-L") interaction
#2]  if we can convert "C=C(C)C.BrBr" to "C([C+](C)C)Br.[Br-]" to test the ("C=C", "Y-L") interaction

class AEtest(unittest.TestCase):
    def setUp(self):
        #(C=C, H-L) 
        self.double = readstring("smi","CC=C.I")
        self.double.addh()
        self.double.connectivity_table = connectivity_table.ConnectivityTable(self.double)
        self.double_sources = segmentation.label_sources(self.double)
        self.double_sinks = segmentation.label_sinks(self.double)
        self.double_products = readstring("smi","C[CH+]C.[I-]")
        self.double_products.addh()
        #(C=C,Y-L) 
        self.double_Y_L= readstring("smi","C=C(C)C.BrBr")
        self.double_Y_L.addh()
        self.double_Y_L.connectivity_table = connectivity_table.ConnectivityTable(self.double_Y_L)
        self.double_Y_L_sources = segmentation.label_sources(self.double_Y_L)
        self.double_Y_L_sinks = segmentation.label_sinks(self.double_Y_L)
        self.double_Y_L_products = readstring("smi","C([C+](C)C)Br.[Br-]")
        self.double_Y_L_products.addh()
    #(C=C, H-L)
    def testdouble(self):
        reaction = reaction_factory.produce_reaction("AE",self.double_sources,self.double_sinks)
        self.assertTrue(reaction.cross_check() == 0.25)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.double,self.double_products)) 
    # (C=C,Y-L)
    def testdouble(self):
        reaction = reaction_factory.produce_reaction("AE",self.double_Y_L_sources,self.double_Y_L_sinks)
        self.assertTrue(reaction.cross_check() == 1)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.double_Y_L,self.double_Y_L_products))

if __name__ == "__main__":
    unittest.main()