#!/usr/bin/python
#coding=utf-8
"""
Structural operations that are more complex than what can be achieved by single OpenBabel function calls.
Usually deals with cross-molecule stuff, which OpenBabel is not generally used to do.
"""

import openbabel
import pybel
import copy
from structure.connectivity_table import ConnectivityTable

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
        new_mol.connectivity_table = ConnectivityTable(new_mol)
    return new_mol

#TODO: Figure out whether we always want order=1 bonds!
def make_bond(start,end,molecule):
    """
    Makes a bond between two atoms by updating connectivity tables in a modified Pybel Molecule object.
    start and end are atom indices.
    
    Since the "combined SMILES" approach we take means that everything is on the 
    SAME Molecule (and thus OBMol) object, we don't need to do strange copying
    shenanigans.
    
    In this program, we can make the assumption that make_bond is always called before break_bond,
    so any strange hypervalent situations (such as H bonded to two atoms) can be left as-is.

    "start" should be a source of electrons, and "end" be a species that is being given electrons,
    otherwise the formal charges won't work out. "start" will have its formal charge increased by 1,
    and end will have its formal charge decreased by 1.
    """
    obmol = molecule.OBMol
    success = obmol.AddBond(start,end,1)
    start_atom = obmol.GetAtom(start)
    end_atom = obmol.GetAtom(end)
    if not success:
        raise ValueError("AddBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
                            %(start,start_atom.GetAtomicNum(),end,end_atom.GetAtomicNum()))
    #TODO: add routine that checks for double bond stuff
    start_atom.SetFormalCharge(start_atom.GetFormalCharge() + 1)
    end_atom.SetFormalCharge(end_atom.GetFormalCharge() - 1)
    #update the connectivity table if the molecule has one - it always should, but callers of this library might not think of that.
    if hasattr(molecule,"connectivity_table"):
        molecule.connectivity_table.add_bond(start,end)

def break_bond(start,end,molecule):
    """
    Removes a bond between two atoms by updating connectivity tables in a modified Pybel Molecule object.
    start and end are atom indices.

    We make the assumption that make_bond is always called before this to avoid the previously-discussed
    problem.

    The "end" atom is considered to be things like hydrogens or leaving groups, which after breaking the bond
    have its formal charge lowered, and the "start" atom would have its formal charge increased.
    """
    obmol = molecule.OBMol
    start_atom = obmol.GetAtom(start)
    end_atom = obmol.GetAtom(end)
    #iterate through all the bonds until we find the one we need to remove
    #this is slow - figure out a better way someday
    found = False
    for bond in openbabel.OBMolBondIter(obmol):
        if (bond.GetBeginAtomIdx() == start and bond.GetEndAtomIdx() == end)  or (bond.GetBeginAtomIdx() == end and bond.GetEndAtomIdx() == start):
            success = obmol.DeleteBond(bond)
            found = True
            if not success:
                raise ValueError("DeleteBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
                                            %(start,start_atom.GetAtomicNum(),end,end_atom.GetAtomicNum()))
            break
    #don't try for/else here. Doesn't work. Don't know why.
    if not found:
        raise ValueError("Bond not found between %s (atomno: %s) and %s (atomno: %s)."
                %(start_atom.GetIdx(),start_atom.GetAtomicNum(),end_atom.GetIdx(),end_atom.GetAtomicNum()))
    start_atom.SetFormalCharge(obmol.GetAtom(start).GetFormalCharge() +1) #TODO: check for double bonds and stuff...
    end_atom.SetFormalCharge(obmol.GetAtom(end).GetFormalCharge() - 1)
    if hasattr(molecule,"connectivity_table"):
        molecule.connectivity_table.remove_bond(start,end)

