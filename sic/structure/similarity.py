#!/usr/bin/python
#coding=utf-8
"""
Handles all operations that have to do with structure similarity.
"""

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
    Finds out whether one state is the same as another by taking the Tanimoto coefficient between the molecules.
    If there is a mapping from state1 to state2 such that the Tanimoto coefficient of each pair is 1,
    then we conclude that state1 is the same as state2.
    """
    #first, return False if lists aren't the same length. Explicit atom balance is still a thing.
    if len(state1) != len(state2):
        return False
    copy_state1 = set(state1)
    copy_state2 = set(state2)
    for st1_mol in copy_state1:
        remove = False
        for st2_mol in copy_state2:
            if tanimoto(st1_mol,st2_mol) == 1:
                print "st2 removed"
                remove = st2_mol
                break
        if remove:
            copy_state2.remove(st2_mol)

    if len(copy_state2) < 1:
        return True
    return False

