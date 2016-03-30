#!/usr/bin/python
#coding=utf-8
"""
Class that identifies a source.
Supersedes the previous ad-hoc data definition, but still in line with specs.
"""

class Source(object):
    """
    Defines a Source. A source has the following data members:

    - subtype : the type of source, e.g. "Y", "C=C"
    - atoms: a dictionary containing all the atoms in the source, with a key
      determining which "piece" of the source the atom corresponds to. This dictionary is mapped to
      the index of each atom rather than the atom object such that we don't have to shift
      references when copying the molecule.
    - molecule: a reference to the molecule object related to this source, for
      convenience in defining structure-modification functions
    """
    def __init__(self,subtype,atoms,molecule):
        """
        Initializes the Source of a particular subtype with a tuple of atom indices
        obtained from the smarts.findall function, and attaches the reference
        to the molecule as well.
        Depending on the subtype, the atoms will be assigned different keys.
        """
        self.subtype = subtype
        self.atoms = {}
        self.molecule = molecule
        for atom in atoms:
            if subtype == "Y":
                self.atoms["Y"] = atom #only one atom total
            elif subtype == "C-":
                self.atoms["C-"] = atom
            #otherwise do some gymnastics
    
    def get_atom(self,key):
        """
        Returns the atom index corresponding to the particular piece of a source. For example, passing
        "Y" as a parameter to this function will return the atom associated with "Y" in a source.
        """
        return self.atoms[key]
