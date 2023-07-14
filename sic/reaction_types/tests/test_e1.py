import unittest
import pka.pka as pka
import reaction_types.reaction_factory as reaction_factory
import structure.similarity as similarity
import structure.struct_ops as struct_ops
import segmentation.segmentation as segmentation
import structure.connectivity_table as connectivity_table
from pybel import readstring

class E1Test(unittest.TestCase):
    """
    Tests whether E1 reactions work by checking whether we can eliminate from a t-butyl chloride that has had the Cl leave.
    """

    def setUp(self):
        self.t_butyl = readstring("smi","C[C+](C)C.[Cl-].CC(C)(C)[O-]")
        self.t_butyl.addh()
        pka.get_all_pka(self.t_butyl)
        self.t_butyl.connectivity_table = connectivity_table.ConnectivityTable(self.t_butyl)
        self.t_butyl_sources = segmentation.label_sources(self.t_butyl)
        self.t_butyl_sinks = segmentation.label_sinks(self.t_butyl)
        self.t_butyl_products = readstring("smi","C=C(C)C.[Cl-].CC(C)(C)O")
        self.t_butyl_products.addh()

    def testCation(self):
        source_to_use = []
        for source in self.t_butyl_sources:
            if source.subtype == "Y" and self.t_butyl.OBMol.GetAtom(source.get_atom("Y")).GetAtomicNum() == 8:
                source_to_use.append(source)
                break
        sink_to_use = []
        for sink in self.t_butyl_sinks:
            if sink.subtype == "C+":
                sink_to_use.append(sink)
                break
        reaction = reaction_factory.produce_reaction("E1",source_to_use,sink_to_use)
        self.assertTrue(reaction.cross_check() == 1.0)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.t_butyl,self.t_butyl_products))

if __name__ == "__main__":
    unittest.main()
