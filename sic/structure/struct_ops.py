#!/usr/bin/python
#coding=utf-8
"""
Structural operations that are more complex than what can be achieved by single OpenBabel function calls.
Usually deals with cross-molecule stuff, which OpenBabel is not generally used to do.
"""

import openbabel
import pybel

#TODO: Figure out whether single_atom flag is enough, and whether we might want to do more per reaction_type. Also figure out how big Ls work here.
def make_bond(start,end,single_atom=False):
	"""
	Makes a bond between two atoms by updating connectivity tables.
        The atom objects that are passed in are the same kind of "atom objects"
        that are used in source/sink identification, i.e. objects that contain
        an "atom" key for the actual Python Atom object, and a "molecule" key
        for the Python Molecule object. Do not change this - the OBMol objects
        cannot be compared.

	If the atoms live on the same molecule, we only need to call AddBond
	on the underlying OBMol object, and all will be well.
	
	However, if the atoms live on two separate molecules, a copy of 
        the source atoms (the "start" atoms) will be made onto the "sink" molecule,
        since by convention the "source" attacks the "sink".
        
        The single_atom flag exists for reactions that only shuffle around a single atom
        rather than making bonds between two molecules. Usually, this means proton transfers,
        but this could be expanded in the future. In this case, we "copy" the sink's atom over
        to the source molecule, since it no longer exists in the sink, and now exists in the source.
        Later, break_bond will delete it from the sink.
	
	In this program, we can make the assumption that make_bond is always called before break_bond,
	so any strange hypervalent situations (such as H bonded to two atoms) can be left as-is.

        If single_atom is false, returns the mapping between indices in the old source molecule
        and the indices of the same atoms in the new composite molecule, such that source-sink data
        can be updated to have the right OBAtom objects.
	"""
        start_mol = start["molecule"]
        end_mol = end["molecule"]
        start_atom = start["atom"]
        end_atom = end["atom"]
	#first check if they are the same object
	if start_mol == end_mol:
		#if they're the same, then it's easy, just make bond
		success = start_mol.OBMol.AddBond(start_atom.idx,end_atom.idx)
		if not success:
			raise ValueError("AddBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
					%(start_atom.idx,start_atom.atomicnum,end_atom.idx,end_atom.atomicnum))
	else:
                if single_atom:
                    #then it gets tough. Copy end atom...
                    end_copy = copy.deepcopy(end_atom.OBAtom)
                    #now put it in to the other molecule
                    start_mol.OBMol.InsertAtom(end_copy) #implicit begin/endModify
                    #and now add the bond
                    success = start_mol.OBMol.AddBond(start_atom.idx,end_copy.OBAtom.GetIdx()) #hopefully this changes...unit tests!
                    if not success:
                            raise ValueError("AddBond failed for bond between %s (atomno: %s) and %s (atomno: %s)."
                                            %(start_atom.idx,start_atom.atomicnum,end_copy.idx,end_copy.atomicnum))
                else:
                    #even tougher. Copy the entirety of the source into the sink, and keep track of there the indices go.
                    #this is done so that we can restore bond connectivity.
                    index_mapping = {}
                    end_mol.OBMol.BeginModify() #so that this doesn't become a huge performance hit
                    for atom in start_mol: #"start" is always source, "end" is always sink.
                        add_atom_success = end_mol.AddAtom(atom.OBAtom)
                        if not add_atom_success:
                            raise ValueError("AddAtom failed to add atom with atomicnum %s, index %s on source."
                                    %(atom.atomicnum,atom.OBAtom.GetIdx()))
                        #atoms are always added at the end, so use len() equivalent 
                        index_mapping[atom.idx] = end_mol.NumAtoms() #indices start by 1, so no -1 here
                    #then now that the mapping is built, add in the bonds...
                    for bond in openbabel.OBMolBondIter(start_mol):
                        add_bond_success = end_mol.AddBond(index_mapping[bond.GetBeginAtomIdx()],index_mapping[bond.GetEndAtomIdx()])
                        #this should work because we copied ALL atoms in start_mol, so we don't need to check
                        if not add_bond_success:
                            raise ValueError("AddBond failed to add bond between %s and %s on source, which are %s and %s on new mol"
                                    %(bond.GetBeginAtomIdx(),bond.GetEndAtomIdx(),index_mapping[bond.GetBeginAtomIdx()],index_mapping[bond.GetEndAtomIdx()])
                    #remember to wipe the molecule from the ReactionState after this.

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


