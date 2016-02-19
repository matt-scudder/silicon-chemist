#!/usr/bin/python
#coding=utf-8
"""
Methods for getting properties of a structure rather than performing operations on it.
Examines the degree of a carbon, for example.
"""

def get_bonds(atom_obj):
    """
    Takes in an atom object (such as source["atoms"]["Y"]), and gets the indices of the atoms it's bound to
    based on the molecule connectivity table.
    Utility function for less typing.
    """
    return atom_obj["molecule"].connectivity_table[atom_obj["atom"].idx]


def get_carbon_degree(s_obj,carbon_label=False):
    """
    Takes in an source/sink object, and gets the number of non-H atoms attached to its carbon.
    Used to determine whether a carbon is primary, secondary, tertiary or methyl.
    If carbon_label is set, uses that to search for the carbon instead of the string "C".
    """
    carb_string = carbon_label if carbon_label else "C"
    carbon_count = 0
    mol = s_obj["atoms"][carb_string]["molecule"]
    has_L = "L" in s_obj["atoms"]
    for bond in get_bonds(s_obj["atoms"][carb_string]):
        if not (mol.OBMol.GetAtom(bond).IsHydrogen()) and (has_L and bond != s_obj["atoms"]["L"]["atom"].idx):
            #H bonds don't count, neither do L if any
            #TODO: Update the above conditional for species other than L that don't count
            carbon_count += 1
    return carbon_count


