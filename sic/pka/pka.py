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
import pybel
from ..utils import get_real_indices

LONE_PAIR_ATOMS = set([7,8,9,15,16,17,35,53]) #atomic numbers that can have lone pairs commonly. Excludes boron and carbon, those are special.

def get_all_pka(molecule):
    """
    For a particular molecule, find the pKa_HA values for all hydrogens,
    and all pKa_BH values for all atoms with lone pairs.
    After running this method, molecule should have its atoms updated with
    a pka_value property, which carries the pKa_HA value for H atoms,
    and pKa_BH value for non-H atoms.
    Certain species (carbocations, multiple bonds) are ignored by this method,
    and require separate checks.
    """
    #the below looks like O(scary), but the lists are small enough that we don't have to care.
    #I know the access method isn't super great, but it's what we lose for the sake of sequential access
    for pka_obj in PKA_CHART:
        pattern = pka_obj.keys()[0] #there's only one key, so this is fine
        smarts = pybel.Smarts(pattern)
        indices = get_real_indices(smarts.findall(molecule))
		for group in indices:
			for atom_idx in group:
				atom = molecule.atoms[atom_idx]
				if atom.atomicnum == 1: #hydrogen
					atom.pka_value = pka_obj[pattern]["pKa_HA"] 
				elif atom.atomicnum in LONE_PAIR_ATOMS:
					atom.pka_value = pka_obj[pattern]["pKa_BH"]
				else:
					#TODO: Update this for C with lone pairs as well as borohydride stuff
					continue
    #return nothing because this just modifies state
