#!/usr/bin/python
#coding=utf-8
"""
Utilities module for dealing with the idiosyncracies of the various libraries
involved in the project.
"""

import copy
import openbabel

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

def make_bond(start,end):
	"""
	Makes a bond between two atoms by updating connectivity tables.
	If the atoms live on the same molecule, we only need to call AddBond
	on the underlying OBMol object, and all will be well.
	
	However, if the atoms live on two separate molecules, a copy of one atom
	(by convention, the "end" atom, which is the sink) will be made and placed in the
	other molecule.
	
	In this program, we can make the assumption that make_bond is always called before break_bond,
	so any strange hypervalent situations (such as H bonded to two atoms) can be left as-is.
	"""
	start_mol = start.OBAtom.GetParent()
	end_mol = end.OBAtom.GetParent()
	#first check if they are the same object
	if start_mol == end_mol:
		#if they're the same, then it's easy, just make bond
		success = start_mol.AddBond(start.OBAtom.GetIdx(),end.OBAtom.GetIdx())
		if not success:
			raise ValueError("AddBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
					%(start.OBAtom.GetIdx(),start.atomicnum,end.OBAtom.GetIdx(),end.atomicnum))
	else:
		#then it gets tough. Copy end atom...
		end_copy = copy.deepcopy(end.OBAtom)
		#now put it in to the other molecule
		start_mol.AddAtom(end_copy)
		#and now add the bond
		success = start_mol.AddBond(start.OBAtom.GetIdx(),end_copy.GetIdx())
		if not success:
			raise ValueError("AddBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
					%(start.OBAtom.GetIdx(),start.atomicnum,end.OBAtom.GetIdx(),end.atomicnum))

def break_bond(start,end,remove_atom):
	"""
	Removes a bond between two atoms by updating connectivity tables.
	The atoms should always live on the same molecule - if they don't, we have problems.
	We make the assumption that make_bond is always called before this to avoid the previously-discussed
	problem.

	Takes in "remove_atom" as a boolean on whether after breaking this bond, the end atom no
	longer has ties to the start molecule. For an example of this, breaking the H-L bond
	means that the H is no longer part of the same molecule as L, but is part of the same molecule
	as whatever just picked off the H, so remove_atom should be passed in as True in this case.
	"""
	#first check if we're on the same molecule
	start_mol = start.OBAtom.GetParent()
	end_mol = end.OBAtom.GetParent()
	if start_mol != end_mol:
		raise ValueError("When breaking bonds, the molecule of the start and end atom should be the same.\
				 Look into why this broke!")
	else:
		if remove_atom:
			#removing an atom removes all the bonds when you call DeleteAtom rather than DestroyAtom
			success = start_mol.DeleteAtom(end.OBAtom)
			if not success:
				raise ValueError("DeleteAtom failed for %s (atomno: %s)."
						%(end.OBAtom.GetIdx(),end.atomicnum))
		else:
			#then things get complicated because bond indexes aren't a thing in obabel yet.
			#iterate through all the bonds until we find the one we need to remove
			#this is slow - figure out a better way someday
			for bond in openbabel.OBMolBondIter(start_mol):
				#TODO: Check if below is sufficient or whether we can't trust that start and end are our start and end
				if bond.GetBeginAtomIdx() == start.OBAtom.GetIdx() and bond.GetEndAtomIdx() == end.OBAtom.GetIdx():
					success = start_mol.DeleteBond(bond)
					raise ValueError("DeleteBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
							%(start.OBAtom.GetIdx(),start.atomicnum,end.OBAtom.GetIdx(),end.atomicnum))
