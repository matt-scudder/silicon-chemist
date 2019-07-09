#!/usr/bin/python
#coding=utf-8
from collections import defaultdict #to make life easier
import openbabel
"""
Contains the ConnectivityTable class,
which carries all the information necessary to determine
a reaction's current state, in order to figure out whether it's closer
to the product.
"""

class ConnectivityTable(object):
    def get_bond_set(self,start_atom,end_atom):
        """
        Gets the frozenset needed as a key for adding or removing
        bonds from the closer_to_product table. Made a function
        since we repeat it so often.
        """
        #first figure out if start_atom or end_atom is an H, this is why we need the mol
        obmol = self.molecule.OBMol
        start_is_H = obmol.GetAtom(start_atom).GetAtomicNum() == 1 
        end_is_H = obmol.GetAtom(end_atom).GetAtomicNum() == 1
        set_to_check = frozenset([start_atom,end_atom])        #base case
        if start_is_H and end_is_H: #special case for NaH      # if the two atoms are H
            set_to_check = frozenset(["H"])
        elif start_is_H or end_is_H:                           # if one is H and the other is not
            if start_is_H:
                set_to_check = frozenset(["H",end_atom])
            else:
                set_to_check = frozenset([start_atom,"H"])
        return set_to_check                                   

    def add_bond(self,start_atom,end_atom,bond_order=1):
        """
        Adds a bond to both connectivity tables.
        """
        #first add to the regular connectivity table
        #since these are defaultdicts of sets, adding twice for double/triple bonds doesn't matter.
        self.connectivity_table[start_atom].add(end_atom)
        self.connectivity_table[end_atom].add(start_atom)
        #now add to the special one for closer to product        
        set_to_add = self.get_bond_set(start_atom,end_atom)
        #TODO: Is this a bug? Should a double bond add two?
        #self.closer_to_product_table[set_to_add] += 1    # Doesn't count multible order donds. It just adds 1 if there is a connection between the two  atoms
        self.closer_to_product_table[set_to_add] += bond_order
        # Since the frozenset is the key, the bonds including hydrogen will have the same key, so the dictionary value will be updated. Ex: {(1,"H"):3}


    def remove_bond(self,start_atom,end_atom):
        """
        Removes a bond from the connectivity table.
        """
        #first check the closer_to_product table, this way we can check bond orders.
        obmol = self.molecule.OBMol
        set_to_remove = self.get_bond_set(start_atom,end_atom)
        num_bonds = self.closer_to_product_table[set_to_remove]
        if num_bonds > 1: #if it's got BO > 1, subtract 1 from closer_to_product table but not from the other
            self.closer_to_product_table[set_to_remove] -= 1
        else: #full remove
            self.connectivity_table[start_atom].remove(end_atom)
            self.connectivity_table[end_atom].remove(start_atom)
            del self.closer_to_product_table[set_to_remove]

    def __init__(self,mol):
        """
        Takes as input a molecule and creates a ConnectivityTable object, which is then attached to
        the molecule taken in as input.

        A ConnectivityTable stores two different tables for different purposes. The first
        is a dict holding the same information as OBMolBondIter, but in a more easily-searchable fashion.
        This way, we can get bond information for several disparate bonds in succession without having to iterate
        through the entire BondIter every time a search is required.
        The format of the dict is as follows:

        {atom_idx : set([bonded_atom_idx,bonded_atom_idx,bonded_atom_idx,bonded_atom_idx])}

        This is reversible, i.e. if atom 1 is bonded to atom 2, then 1 will be present with 2 among its bonds,
        and 2 will be present with 1 among its bonds.

        A secondary connectivity table is also created for the purpose of computing the closer_to_product function,
        which holds frozensets of bonds as keys, and the order of the bond as values. This is done because an OBMol object
        stores a single bond of multiple order, if there is a double or triple bond. Also, for the sake of closer_to_product,
        we don't really care about hydrogen atom indices, just that they are hydrogens. Therefore, we spec the second connectivity
        table to be as follows:

        {frozenset([1,"H"]):3,frozenset([2,3]):1}

        Where the first key implies that there are 3 bonds from atom 1 to H atoms, and the second key implies a single bond from atom 2 to atom 3.
        Since bonds are bidirectional, we store a set.
        """
        self.molecule = mol #pointer for making add_bond and remove_bond easier
        self.connectivity_table = defaultdict(set) #defaultdicts are the best
        self.closer_to_product_table = defaultdict(int) #defaults to 0
        for bond in openbabel.OBMolBondIter(self.molecule.OBMol): #TODO: Does this count a double-bond twice? (See TODO in add_bond()) 
            start_atom = bond.GetBeginAtomIdx()                   #### "openbabel.OBMolBondIter" counts the "openbabel.OBBond". If its a double bond, it stores this info in the "openbabel.OBBond" object 
            end_atom = bond.GetEndAtomIdx()                       #### we can get the bond order using the "bond.GetBondOrder()" fnction
            bond_order = bond.GetBondOrder()
            self.add_bond(start_atom,end_atom,bond_order)

    def __repr__(self):
        """
        Returns a string representation of this object for easy debugging.
        """
        return "ConnectivityTable<Table:{}>".format(self.connectivity_table)

    def get_atoms_bonded(self,atom):
        """
        Gets all the atoms bonded to a particular atom. Takes as input an atom index, and returns
        a set of all atoms bonded to a particular atom.
        """
        return self.connectivity_table[atom]
    
    def bond_exists(self,start_atom,end_atom):
        """
        Takes as input a start atom and an end atom and checks whether a bond between them exists.
        """
        return end_atom in self.connectivity_table[start_atom] and start_atom in self.connectivity_table[end_atom] #both sides should alwaays be true

    def get_bond_degree(self,start_atom,end_atom):
        """
        Takes as input a start atom and an end atom and checks the degree of the bond.
        If the bond is not present, returns 0.
        """
        set_to_check = self.get_bond_set(start_atom,end_atom)
        if set_to_check not in self.closer_to_product_table:
            return 0
        else:
            return self.closer_to_product_table[set_to_check]

    def get_all_unidirectional_bonds(self):
        """
        Returns all the unidirectional bonds, i.e. the bonds as listed in the closer_to_product_table.
        """
        return self.closer_to_product_table.keys()

    def get_bond_degree(self,bondset):
        """
        Takes a frozenset in the format of the closer_to_product table, and gets the degree of the
        according bond. If one of the members is H, returns the number of bonds to hydrogen atoms
        originating at the other atom.
        """
        return self.closer_to_product_table[bondset]


    def get_bonds_to_atomicnum(self,atom_idx,atomicnum):
        """
        Base function that gets the number of bonds from an atom to a particular atomic number.
        """
        bond_count = 0
        for bonded_atom_idx in self.table[atom_idx]:
            bonded_atom_obj = self.molecule.GetAtom(bonded_atom_idx)
            if bonded_atom_obj.GetAtomicNum() == atomicnum:
                bond_count += 1
        return bond_count

    def get_hydrogens(self,atom_idx):
        """
        Returns the number of hydrogens attached to a particular atom.
        Convenience function.
        """
        return get_bonds_to_atomicnum(atom_idx,1)
