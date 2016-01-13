#!/usr/bin/python
#coding=utf-8
"""
Structural operations that are more complex than what can be achieved by single OpenBabel function calls.
Usually deals with cross-molecule stuff, which OpenBabel is not generally used to do.
"""

import openbabel

#TODO: Figure out what to do when you "bond" two molecules together by one link...
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

def break_bond(start,end):
	"""
	Removes a bond between two atoms by updating connectivity tables.
	The atoms should always live on the same molecule - if they don't, we have problems.
	We make the assumption that make_bond is always called before this to avoid the previously-discussed
	problem.
	The "end" atom is considered to be things like hydrogens or leaving groups, which after breaking the bond
	no longer have any ties to the molecule and may be deleted because of it.
	"""
	#first check if we're on the same molecule
	start_mol = start.OBAtom.GetParent()
	end_mol = end.OBAtom.GetParent()
	if start_mol != end_mol:
		raise ValueError("When breaking bonds, the molecule of the start and end atom should be the same.\
				 Look into why this broke!")
	else:
		#iterate through all the bonds until we find the one we need to remove
		#this is slow - figure out a better way someday
		for bond in openbabel.OBMolBondIter(start_mol):
			#TODO: Check if below is sufficient or whether we can't trust that start and end are our start and end
			if bond.GetBeginAtomIdx() == start.idx and bond.GetEndAtomIdx() == end.idx:
				success = start_mol.DeleteBond(bond)
				if not success:
					raise ValueError("DeleteBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
							%(start.OBAtom.GetIdx(),start.atomicnum,end.OBAtom.GetIdx(),end.atomicnum))
		if end.implicitvalence < 1:
			#if no bonds left, delete the atom
			#this may be lies, so make sure this actually works
			success = start_mol.DeleteAtom(end.OBAtom)
			if not success:
				raise ValueError("DeleteAtom failed for %s (atomno: %s)."
						%(end.OBAtom.GetIdx(),end.atomicnum))


