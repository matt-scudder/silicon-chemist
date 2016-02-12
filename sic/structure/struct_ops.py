#!/usr/bin/python
#coding=utf-8
"""
Structural operations that are more complex than what can be achieved by single OpenBabel function calls.
Usually deals with cross-molecule stuff, which OpenBabel is not generally used to do.
"""

import openbabel
import pybel

def add_bond_connectivity_table(start_atom,end_atom,table):
    """
    Adds a bond to the connectivity table.
    This is called several times throughout struct_ops, so it's a function.
    """
    if start_atom not in table:
        table[start_atom] = set([end_atom])
    else:
        table[start_atom].add(end_atom)
    if end_atom not in table:
        table[end_atom] = set([start_atom])
    else:
        table[end_atom].add(start_atom)

def remove_bond_connectivity_table(start_atom,end_atom,table):
    """
    Removes a bond from the connectivity table.
    If attempting to remove a bond that doesn't exist, does not change
    the connectivity table.
    This is called several times throughout struct_ops, so it's a function.
    """
    #find start atom, remove end_atom from its set
    #if start atom is not in the table, no worries - it already "doesn't have any bonds",
    #so removing them isn't relevant and we can fail silently.
    if start_atom in table:
        table[start_atom].remove(end_atom)
    if end_atom in table: #end_atom should be present iff start_atom is present, but might as well clean up.
        table[end_atom].remove(start_atom)

def generate_connectivity_table(mol):
    """
    Takes as input a molecule and attaches a connectivity table to said molecule,
    which is a dict holding the same information as OBMolBondIter, but in a more easily-searchable fashion.
    This way, we can get bond information for several disparate bonds in succession without having to iterate
    through the entire BondIter every time a search is required.
    The format of the dict is as follows:

    {atom_idx : set([bonded_atom_idx,bonded_atom_idx,bonded_atom_idx,bonded_atom_idx])}

    This is reversible, i.e. if atom 1 is bonded to atom 2, then 1 will be present with 2 among its bonds,
    and 2 will be present with 1 among its bonds.
    """
    table = {}
    for bond in openbabel.OBMolBondIter(mol.OBMol):
        start_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        add_bond_connectivity_table(start_atom,end_atom,table)
    mol.connectivity_table = table


def copy_molecule(mol):
    """
    Copies a molecule by using the operator= of OBMol as well as inserting references
    to the objects added to Molecule objects by our program in the new Molecule object.
    Used such that atom indexing for "closer to product" mappings can be kept
    consistent throughout a mechanism.
    """
    intermediate = openbabel.OBMol(mol.OBMol) #Molecule's constructor only takes OBMol objects, so first copy the one from the original
    #using the OBMol constructor, which copies all atoms
    new_mol = pybel.Molecule(intermediate) #add it to a new Molecule object
    #copy the properties only in the case when they aren't there, to prevent strange bugs
    #when people call this function as a utility for other purposes.
    if hasattr(mol,"pka_index"):
        new_mol.pka_index = mol.pka_index
    if hasattr(mol,"connectivity_table"):
        new_mol.connectivity_table = mol.connectivity_table
    return new_mol

#TODO: Figure out whether we always want order=1 bonds!
def make_bond(start,end):
    """
    Makes a bond between two atoms by updating connectivity tables.
    The atom objects that are passed in are the same kind of "atom objects"
    that are used in source/sink identification, i.e. objects that contain
    an "atom" key for the actual Python Atom object, and a "molecule" key
    for the Python Molecule object. Do not change this - the OBMol objects
    cannot be compared.
    
    Since the "combined SMILES" approach we take means that everything is on the 
    SAME Molecule (and thus OBMol) object, we don't need to do strange copying
    shenanigans.
    
    In this program, we can make the assumption that make_bond is always called before break_bond,
    so any strange hypervalent situations (such as H bonded to two atoms) can be left as-is.

    This function does not use Atom objects in order to make shifting atom references
    less necessary.
    """
    #because of infrastructure changes, start_mol == end_mol always
    start_mol = start["molecule"]
    end_mol = end["molecule"]
    start_atom = start["atom"]
    end_atom = end["atom"]
    success = start_mol.OBMol.AddBond(start_atom.idx,end_atom.idx,1)
    #TODO: add routine that checks for double bond stuff
    start_atom.OBAtom.SetFormalCharge(start_mol.OBMol.GetAtom(start_atom.idx).GetFormalCharge() + 1)
    if end_atom.atomicnum != 1: #hydrogen behaves oddly w.r.t. formal charges
        end_atom.OBAtom.SetFormalCharge(end_mol.OBMol.GetAtom(end_atom.idx).GetFormalCharge() - 1)
    #update the connectivity table if the molecule has one - it always should, but callers of this library might not think of that.
    if not success:
        raise ValueError("AddBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
                            %(start_atom.idx,start_atom.atomicnum,end_atom.idx,end_atom.atomicnum))
    if hasattr(start_mol,"connectivity_table"):
        add_bond_connectivity_table(start_atom.idx,end_atom.idx,start_mol) 

def break_bond(start,end):
    """
    Removes a bond between two atoms by updating connectivity tables.
    The atom objects that are passed in are the same kind of "atom objects"
    that are used in source/sink identification, i.e. objects that contain
    an "atom" key for the actual Python Atom object, and a "molecule" key
    for the Python Molecule object. Do not change this - the OBMol objects
    cannot be compared.
    The atoms should always live on the same molecule - if they don't, we have problems.
    We make the assumption that make_bond is always called before this to avoid the previously-discussed
    problem.
    The "end" atom is considered to be things like hydrogens or leaving groups, which after breaking the bond
    may no longer have any ties to the molecule and may be deleted because of it.
    """
    #first check if we're on the same molecule
    start_mol = start["molecule"]
    end_mol = end["molecule"]
    start_atom = start["atom"]
    end_atom = end["atom"]
    #iterate through all the bonds until we find the one we need to remove
    #this is slow - figure out a better way someday
    found = False
    for bond in openbabel.OBMolBondIter(start_mol.OBMol):
        if (bond.GetBeginAtomIdx() == start_atom.idx and bond.GetEndAtomIdx() == end_atom.idx)  or (bond.GetBeginAtomIdx() == end_atom.idx and bond.GetEndAtomIdx() == start_atom.idx):
            success = start_mol.OBMol.DeleteBond(bond)
            found = True
            if not success:
                raise ValueError("DeleteBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
                                            %(start_atom.OBAtom.GetIdx(),start_atom.atomicnum,end_atom.OBAtom.GetIdx(),end_atom.atomicnum))
            break
    #don't try for/else here. Doesn't work. Don't know why.
    if not found:
        raise ValueError("Bond not found between %s (atomno: %s) and %s (atomno: %s)."
                %(start_atom.OBAtom.GetIdx(),start_atom.atomicnum,end_atom.OBAtom.GetIdx(),end_atom.atomicnum))
    start_atom.OBAtom.SetFormalCharge(start_mol.OBMol.GetAtom(start_atom.idx).GetFormalCharge() -1) #TODO: check for double bonds and stuff...
    if hasattr(start_mol,"connectivity_table"):
        remove_bond_connectivity_table(start_atom.idx,end_atom.idx,start_mol)

