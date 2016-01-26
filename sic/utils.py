#!/usr/bin/python
#coding=utf-8
"""
Utilities module for dealing with the idiosyncracies of the various libraries
involved in the project.
"""

import copy

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
    return state.write("smiles").rstrip()
