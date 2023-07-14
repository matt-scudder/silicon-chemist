"""
Class that identifies a sink.
Supersedes the previous ad-hoc data definition, but still in line with specs.
"""

class Sink(object):
    """
    Defines a Sink. A sink has the following data members:

    - subtype : the type of sink, e.g. "H-L", "C-L"
    - atoms: a dictionary containing all the atoms in the sink, with a key
      determining which "piece" of the sink the atom corresponds to. This dictionary is mapped to
      the index of each atom rather than the atom object such that we don't have to shift
      references when copying the molecule.
    - molecule: a reference to the molecule object related to this sink, for
      convenience in defining structure-modification functions
    """
    def __init__(self,subtype,atoms,molecule):
        """
        Initializes the Sink of a particular subtype with a tuple of atom indices
        obtained from the smarts.findall function, and attaches the reference
        to the molecule as well.
        Depending on the subtype, the atoms will be assigned different keys.
        """
        self.subtype = subtype
        self.atoms = {}
        self.molecule = molecule
        for atom in atoms:
            if subtype == "H-L":
                if molecule.OBMol.GetAtom(atom).GetAtomicNum() == 1: #to avoid molecule.atoms which is a list comprehension on OBMol
                    self.atoms["H"] = atom
                else:
                    self.atoms["L"] = atom #H-H we pretend not to care about
            elif subtype == "C+":
                self.atoms["C+"] = atom #only one atom to care about
            elif subtype == "C-L":
                #TODO: Figure out what to do if L is also a C, though C-C bonds don't break often.
                if molecule.OBMol.GetAtom(atom).GetAtomicNum() == 6:
                    self.atoms["C"] = atom
                else:
                    self.atoms["L"] = atom
            elif subtype == "ewg(O)-C=C" or subtype == "ewg(N)-C=C": # Each "ewg-C=C" source has 4 atoms:(F, S, C1 , C2). F = first, S= second , C1 and C2 are the C-atoms of the double bond
                if molecule.OBMol.GetAtom(atoms[1]).GetAtomicNum() == 8:
                    self.atoms["F"] = atoms[1]
                    self.atoms["S"] = atom
                else:
                    self.atoms["F"] = atom
                    self.atoms["S"] = atoms[1]
                
                self.atoms["C1"] = atoms[2]
                self.atoms["C2"] = atoms[3]
                break
            elif subtype == "Z=C":
                if molecule.OBMol.GetAtom(atom).GetAtomicNum() == 6:
                    self.atoms["C"] = atom
                else:
                    self.atoms["Z"] = atom
            elif subtype == "Y-L":
                self.atoms["Y"] = atom
                self.atoms["L"] = atoms[1]
                break
        


            #otherwise do some gymnastics
    
    def get_atom(self,key):
        """
        Returns the atom index corresponding to the particular piece of a sink. For example, passing
        "C" as a parameter to this function will return the atom associated with "C" in a sink.
        """
        return self.atoms[key]

    def __repr__(self):
        """
        Returns a string representation of this object for easy debugging.
        """
        return "Sink<Type:{},Atoms:{}>".format(self.subtype,self.atoms)
