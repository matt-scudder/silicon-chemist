#!/usr/bin/python
#coding=utf-8
"""
This part of the program handles pKa determination.
It stores the pKa chart as a list of SMARTS strings, 
with each entry having a pKa_HA and a pKa_bH.
To identify the pKa of all atoms in a molecule for
which it makes sense to do so (hydrogens for pKa_HA
and O,S,N, and halogens for pKa_BH), this module
will go through the species in the pKa
chart and attempt to map them to the parts of the molecule.
"""
from pka_chart import PKA_CHART
from openbabel import pybel

LONE_PAIR_ATOMS = set([6,7,8,9,15,16,17,35,53]) #atomic numbers that can have lone pairs commonly. Excludes boron , this is special.

def get_all_pka(molecule):
    """
    For a particular molecule, find the pKa_HA values for all hydrogens,
    and all pKa_BH values for all atoms with lone pairs.
    After running this method, molecule should have a pka_index attribute,
    which maps atom indices to their pKa values, specifically
    the pKa_HA value for H atoms, and pKa_BH value for non-H atoms.
    Certain species (carbocations, multiple bonds) are ignored by this method,
    and require separate checks.
    """
    #the below looks like O(scary), but the lists are small enough that we don't have to care.
    #I know the access method isn't super great, but it's what we lose for the sake of sequential access
    molecule.pka_index = {} #This is cleared and recalculated at each time in order to keep it consistent with changes in bonds
    for pka_obj in PKA_CHART:
        pattern = pka_obj.keys()[0] #there's only one key, so this is fine
        smarts = pybel.Smarts(pattern)
        indices = smarts.findall(molecule)
        obmol = molecule.OBMol
        for group in indices:
            h_atom = False
            for atom_idx in group:
                atom = obmol.GetAtom(atom_idx)
                if atom.GetAtomicNum() == 1: # get all hydrogens, they should all be bonded to the same atom.
                    h_atom = atom_idx
                    molecule.pka_index[atom_idx] = pka_obj[pattern]["pKa_HA"] 
            if h_atom:
                #get atom bonded to H if there is one
                bonded_to_h = molecule.connectivity_table.get_atoms_bonded(h_atom) 
                molecule.pka_index[next(iter(bonded_to_h))] = pka_obj[pattern]["pKa_BH"] #this should ONLY have one element in the set, so it works out
            else:
                for atom_idx in group:
                    if atom.GetAtomicNum() in LONE_PAIR_ATOMS:
                        molecule.pka_index[atom_idx] = pka_obj[pattern]["pKa_BH"]
    #return nothing because this just modifies state

def get_pka(atom,molecule):
    """
    Takes in an atom index and its corresponding Molecule,
    and returns the pKa value associated to the atom by checking the Molecule's pKa index.

    The Molecule object is required because OBMol (SWIG proxy object) does not seem to have persistent attributes,
    i.e. they get cleared every time the OBMol object is accessed anew, likely because of descriptor use.

    If the pKa value is either None or not present, None will be returned (since the result is equivalent).

    We assume that pka_index is an attribute of molecule. If this is not the case, then this method will
    raise an AttributeError.
    """
    index = atom
    if molecule.pka_index.has_key(index):
        return molecule.pka_index[index]
    else:
        return None
