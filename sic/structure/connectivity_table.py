#!/usr/bin/python
#coding=utf-8
from collections import defaultdict #to make life easier
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
        set_to_check = frozenset([start_atom,end_atom]) #base case
        if start_is_H and end_is_H: #special case for NaH
            set_to_check = frozenset(["H"])
        elif start_is_H or end_is_H:
            if start_is_H:
                set_to_check = frozenset(["H",end_atom])
            else:
                set_to_check = frozenset([start_atom,"H"])
        return set_to_check

    def add_bond(self,start_atom,end_atom):
        """
        Adds a bond to both connectivity tables.
        """
        #first add to the regular connectivity table
        self.connectivity_table[start_atom].add(end_atom)
        self.connectivity_table[end_atom].add(start_atom)
        #now add to the special one for closer to product
        set_to_add = self.get_bond_set(start_atom,end_atom)
        self.closer_to_product_table[set_to_add] += 1


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

        Upon initialization, a ConnectivityTable attaches itself to the input molecule, and sets a reference within itself to the molecule.
        """
        self.molecule = mol #pointer for making add_bond and remove_bond easier
        self.connectivity_table = defaultdict(set) #defaultdicts are the best
        self.closer_to_product_table = defaultdict(int) #defaults to 0
        for bond in openbabel.OBMolBondIter(self.molecule.OBMol):
            start_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            self.add_bond(start_atom,end_atom,self.molecule)
        mol.connectivity_table = self
