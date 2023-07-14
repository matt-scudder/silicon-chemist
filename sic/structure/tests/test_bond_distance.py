"""
Tests whether bond distance is 0 for molecules that are equal,
and that it goes down when you get closer to a product.
"""
from structure import struct_ops, properties, connectivity_table
from segmentation import segmentation
import unittest
from openbabel.pybel import readstring

class GroupAttachmentTest(unittest.TestCase):
    def setUp(self):
        t_butyl_reactants = readstring("smi","CC(C)(C)O.Cl")
        t_butyl_products = readstring("smi","CC(C)(C)Cl.O")
        t_butyl_reactants.addh()
        t_butyl_products.addh()
        t_butyl_reactants.connectivity_table = connectivity_table.ConnectivityTable(t_butyl_reactants)
        t_butyl_products.connectivity_table = connectivity_table.ConnectivityTable(t_butyl_products)
        print(t_butyl_reactants.connectivity_table.connectivity_table)
        #generate t_butyl_mapping here because we're testing the function... 
        self.t_butyl_mapping = properties.get_mapping(t_butyl_reactants,t_butyl_products)
        self.t_butyl_reactants = t_butyl_reactants
        self.t_butyl_products = t_butyl_products
        self.t_butyl_sources = segmentation.label_sources(t_butyl_reactants)
        self.t_butyl_sinks = segmentation.label_sinks(t_butyl_reactants)
        #now same thing for nitriles, which are a problem
        nitrile_reactants= readstring("smi","CC#N.C[O-]")
        nitrile_products = readstring("smi","[CH2-]C#N.CO")
        nitrile_reactants.removeh()
        nitrile_products.removeh()
        nitrile_reactants.addh()
        nitrile_products.addh()
        nitrile_reactants.connectivity_table = connectivity_table.ConnectivityTable(nitrile_reactants)
        nitrile_products.connectivity_table = connectivity_table.ConnectivityTable(nitrile_products)
        #generate nitrile_mapping here because we're testing the function... 
        self.nitrile_mapping = properties.get_mapping(nitrile_reactants,nitrile_products)
        self.nitrile_reactants = nitrile_reactants
        self.nitrile_products = nitrile_products
        self.nitrile_sources = segmentation.label_sources(nitrile_reactants)
        self.nitrile_sinks = segmentation.label_sinks(nitrile_reactants)


    def testAttachment(self):
        #do proton_transfer between OH and HCl
        #NOTE: the numbers below are empirical based on how segmentation operates...
        Y = 5
        H = 17 
        L = 6
        print("Y is idx %s, L is idx %s" % (Y,L))
        first_distance = properties.get_bond_distance(self.t_butyl_reactants,self.t_butyl_products,self.t_butyl_mapping)
        print("Initial distance: %s" % first_distance)
        print(self.t_butyl_reactants.write("can"))
        print("connectivity table crap: {}".format(self.t_butyl_reactants.connectivity_table.get_bond_set(Y,H)))
        struct_ops.make_bond(Y,H,self.t_butyl_reactants) #self.t_butyl_reactants is a Molecule object
        print(self.t_butyl_reactants.write("can"))
        struct_ops.break_bond(H,L,self.t_butyl_reactants) #self.t_butyl_reactants is a Molecule object
        print(self.t_butyl_reactants.write("can"))
        second_distance = properties.get_bond_distance(self.t_butyl_reactants,self.t_butyl_products,self.t_butyl_mapping)
        print("After bond-making and bond-breaking: %s" % second_distance)
        self.assertTrue(second_distance < first_distance)

    def testNitriles(self):
        """
        Tests nitriles, which somehow break the bond distance system.
        """
        Y = 5 
        H = 6
        L = 1
        print("Atoms in reactant:")
        for atom in self.nitrile_reactants.atoms:
            print("Atom: idx: {}, atomicnum: {}".format(atom.idx,atom.atomicnum))
        print("Atoms in reactant:")
        for atom in self.nitrile_products.atoms:
            print("Atom: idx: {}, atomicnum: {}".format(atom.idx,atom.atomicnum))

        first_distance = properties.get_bond_distance(self.nitrile_reactants,self.nitrile_products,self.nitrile_mapping)
        print("Initial distance: %s" % first_distance)
        print(self.nitrile_reactants.write("can"))
        print("connectivity table crap: {}".format(self.nitrile_reactants.connectivity_table.get_bond_set(Y,H)))
        struct_ops.make_bond(Y,H,self.nitrile_reactants) #self.nitrile_reactants is a Molecule object
        print(self.nitrile_reactants.write("can"))
        struct_ops.break_bond(H,L,self.nitrile_reactants) #self.nitrile_reactants is a Molecule object
        print(self.nitrile_reactants.write("can"))
        second_distance = properties.get_bond_distance(self.nitrile_reactants,self.nitrile_products,self.nitrile_mapping)
        print("After bond-making and bond-breaking: %s" % second_distance)
        self.assertTrue(second_distance < 1)


if __name__ == "__main__":
    unittest.main()

