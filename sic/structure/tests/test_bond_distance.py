#!/usr/bin/python
#coding=utf-8
"""
Tests whether bond distance is 0 for molecules that are equal,
and that it goes down when you get closer to a product.
"""
import structure.struct_ops as struct_ops
import structure.similarity as similarity
import structure.properties as properties
import segmentation.segmentation as segmentation
import structure.connectivity_table as connectivity_table
import unittest
from pybel import readstring

class GroupAttachmentTest(unittest.TestCase):
    def setUp(self):
        reactants = readstring("smi","CC(C)(C)O.Cl")
        products = readstring("smi","CC(C)(C)Cl.O")
        reactants.addh()
        products.addh()
        reactants.connectivity_table = connectivity_table.ConnectivityTable(reactants)
        products.connectivity_table = connectivity_table.ConnectivityTable(products)
        #generate mapping here because we're testing the function... 
        self.mapping = properties.get_mapping(reactants,products)
        self.reactants = reactants
        self.products = products
        self.sources = segmentation.label_sources(reactants)
        self.sinks = segmentation.label_sinks(reactants)

    def testAttachment(self):
        #do proton_transfer between OH and HCl
        #NOTE: the numbers below are empirical based on how segmentation operates...
        Y = self.sources[0].get_atom("Y")
        H = self.sinks[1].get_atom("H")
        L = self.sinks[1].get_atom("L")
        print "Y is idx %s, L is idx %s" % (Y,L)
        first_distance = properties.get_bond_distance(self.reactants,self.products,self.mapping)
        print "Initial distance: %s" % first_distance
        print self.reactants.write("can")
        print "connectivity table crap: {}".format(self.reactants.connectivity_table.get_bond_set(Y,H))
        struct_ops.make_bond(Y,H,self.reactants) #self.reactants is a Molecule object
        print self.reactants.write("can")
        struct_ops.break_bond(H,L,self.reactants) #self.reactants is a Molecule object
        print self.reactants.write("can")
        second_distance = properties.get_bond_distance(self.reactants,self.products,self.mapping)
        print "After bond-making and bond-breaking: %s" % second_distance
        self.assertTrue(second_distance < first_distance)

if __name__ == "__main__":
    unittest.main()

