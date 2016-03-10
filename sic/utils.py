#!/usr/bin/python
#coding=utf-8
"""
Utilities module for dealing with the idiosyncracies of the various libraries
involved in the project.
"""
from pybel import Molecule, Atom
from segmentation.sink import Sink
from segmentation.source import Source
from openbabel import OBMolBondIter, OBAtom

def get_real_indices(indices):
    """
    OpenBabel's smarts module returns atom indices starting from 1, not 0
    (see http://openbabel.org/docs/2.3.1/UseTheLibrary/Python_PybelAPI.html ),
    and this appears to be a design choice:
    http://openbabel.org/wiki/Atom_Indexing
    However, the only way to access an atom in a Molecule object is by accessing
    mol.atoms, which is obviously 0-indexed. This is irritating and kind of
    pointless, so this function aims to take in the garbage that we can't use
    without fixing the indices, and changes all indices to something we
    can use directly. That way we avoid future bugs when someone forgets
    to take index-1 instead of index.

    NOTE: While smarts.findall() returns a list of tuples, this function
    will return a list of lists, because the overhead of creating these
    short lists for the size of molecules we're using isn't significant enough
    to do this "right" and keeping type.
    """
    #first return lists that have index-1, then copy them to a new list using list-processing
    return [map(lambda x: x-1,x) for x in indices]

def write_mol(state):
    """
    Pybel for some reason thinks that appending \t\n to every string
    it prints out is fun.
    I don't.
    Therefore this function.
    """
    return state.write("can").rstrip()

def deepcopy_ignoring_mol(item,new_mol):
    """
    Copies our source or sink objects while ignoring Molecule objects,
    since they cannot be deepcopied due to the fact that SWIG gdoes
    not support it directly.
    """
    new_mol_atoms = new_mol.atoms #for copying Atom objects
    if isinstance(item,dict):
        return dict(map(lambda kv: (kv[0], deepcopy_ignoring_mol(kv[1],new_mol)), item.items()))
    elif isinstance(item,Molecule):
        return new_mol
    elif isinstance(item,Atom):
        return new_mol_atoms[item.idx-1]
    elif isinstance(item,Source):
        item.molecule = new_mol
        return item
    elif isinstance(item,Sink):
        item.molecule = new_mol
        return item
    else:
        return item

def shift_molecule_references(class_list,new_mol):
    """
    Shifts the molecule references of lists of either source or sink objects
    to the molecule in new_mol by creating a new list, creating copies
    of the source/sink objects, and returns the new list.
    """
    new_list = []
    for generic_class in class_list:
        new_list.append(deepcopy_ignoring_mol(generic_class,new_mol))
    return new_list 

def write_all_bonds(mol):
    """
    Returns a string of all bonds in a Molecule object with the following format:
    
    (start_idx,end_idx;order)

    Useful for debugging.
    """
    all_bonds = ""
    for bond in OBMolBondIter(mol.OBMol):
        all_bonds += ("(%s,%s;%s)\n"%(bond.GetBeginAtomIdx(),bond.GetEndAtomIdx(),bond.GetBO()))
    return all_bonds



