#!/usr/bin/python
#coding=utf-8
"""
Structural operations that are more complex than what can be achieved by single OpenBabel function calls.
Usually deals with cross-molecule stuff, which OpenBabel is not generally used to do.
"""

import openbabel
import pybel

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
    start_mol = start["molecule"]
    end_mol = end["molecule"]
    start_atom = start["atom"]
    end_atom = end["atom"]
    success = start_mol.OBMol.AddBond(start_atom.idx,end_atom.idx,1)
    #TODO: add routine that checks for double bond stuff
    start_atom.OBAtom.SetFormalCharge(start_mol.OBMol.GetAtom(start_atom.idx).GetFormalCharge() + 1)
    if end_atom.atomicnum != 1: #hydrogen behaves oddly w.r.t. formal charges
        end_atom.OBAtom.SetFormalCharge(end_mol.OBMol.GetAtom(end_atom.idx).GetFormalCharge() - 1)
    if not success:
        raise ValueError("AddBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
                            %(start_atom.idx,start_atom.atomicnum,end_atom.idx,end_atom.atomicnum))

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
    if end_atom.valence < 1:
        #if no bonds left, delete the atom
        #this may be lies, so make sure this actually works
        success = start_mol.OBMol.DeleteAtom(end_atom.OBAtom)
        if not success:
            raise ValueError("DeleteAtom failed for %s (atomno: %s)."
                            %(end_atom.OBAtom.GetIdx(),end_atom.atomicnum))
