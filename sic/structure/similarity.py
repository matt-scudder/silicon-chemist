"""
Handles all operations that have to do with structure similarity.
"""

from sic.utils import write_mol


def tanimoto(mol1,mol2):
    """
    Returns the tanimoto coefficient between two molecules, using the FP3 fingerprint.
    Because FP2 ignores single-atom C,N,O fragments, and these are often of crucial importance
    to determining whether we're at product or not, FP2 is not used. Other fingerprints may be
    used if FP3 is insufficient.
    
    If the coefficient is 1, then the two molecules are definitely the same.
    """
    return mol1.calcfp(fptype="fp3") | mol2.calcfp(fptype="fp3")

def is_same_state(state1,state2):
    """
    Returns True if both states are composed of the same SMILES strings.
    If we can't tell the molecules apart based on SMILES, we can conclude that they are the same -
    since SMILES is our input, that is the highest level of resolution that we get.
    """
    state1_strings = set(normalize_mols(state1))
    state2_strings = set(normalize_mols(state2))
    return state1_strings == state2_strings

def is_same_molecule(mol1,mol2):
    """
    Determines whether two Molecule objects are the same thing.
    """
    return write_mol(mol1) == write_mol(mol2)

def normalize_mols(mols):
    """
    Normalizes a list of Molecules in "can" format.
    """
    return [write_mol(x) for x in mols]
