#!/usr/bin/python
#coding=utf-8
import unittest
import pka.pka as pka
import reaction_types.reaction_factory as reaction_factory
import structure.similarity as similarity
import structure.struct_ops as struct_ops
import segmentation.segmentation as segmentation
import structure.connectivity_table as connectivity_table
from pybel import readstring

class EBTest(unittest.TestCase):
    """
    Tests whether EB reactions work by checking whether we can eliminate from 3-chloro-propanenitrile that has been deprotonated.
    """

    def setUp(self):
        self.nitrile = readstring("smi","N#C[CH-]CCl")
        self.nitrile.addh()
        pka.get_all_pka(self.nitrile)
        self.nitrile.connectivity_table = connectivity_table.ConnectivityTable(self.nitrile)
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
        self.assertTrue(reaction.cross_check() == 1.0)
        reaction.rearrange()
        self.assertTrue(similarity.is_same_molecule(self.nitrile,self.nitrile_products))

if __name__ == "__main__":
    unittest.main()
