import unittest

from openbabel.pybel import readstring

from sic.pka import pka
from sic.reaction_types import reaction_factory
from sic.segmentation import segmentation
from sic.structure import similarity, connectivity_table

class EBTest(unittest.TestCase):
    """
    Tests whether EB reactions work by checking whether we can eliminate from 3-chloro-propanenitrile that has been deprotonated.
    """

    def setUp(self):
        self.nitrile = readstring("smi","N#C[CH-]CCl")
        self.nitrile.addh()
        self.nitrile.connectivity_table = connectivity_table.ConnectivityTable(self.nitrile)
        pka.get_all_pka(self.nitrile)
        self.nitrile_sources = segmentation.label_sources(self.nitrile)
        self.nitrile_sinks = segmentation.label_sinks(self.nitrile)
        self.nitrile_products = readstring("smi","C=CC#N.[Cl-]")
        self.nitrile_products.addh()

    def testCation(self):
        source_to_use = []
        for source in self.nitrile_sources:
            if source.subtype == "C-":
                source_to_use.append(source)
                break
        sink_to_use = []
        for sink in self.nitrile_sinks:
            if sink.subtype == "C-L":
                sink_to_use.append(sink)
                break
        reaction = reaction_factory.produce_reaction("EB",source_to_use,sink_to_use)
        self.assertEquals(reaction.cross_check(), 1.0)
        reaction.rearrange()
        self.assertEquals(*similarity.normalize_mols([self.nitrile,self.nitrile_products]))

if __name__ == "__main__":
    unittest.main()
