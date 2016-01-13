#!/usr/bin/python
#coding=utf-8
"""
While the sources and sinks "modules" contain the actual expressions used
for matching, this module contains the logic used to actually perform the
matching on a molecule and generate the data structures other parts of SiC³
will use to make sure things are rearranged properly when a reaction occurs.
"""

from sources import SOURCES
from sinks import SINKS
import pybel
from ..utils import get_real_indices 

#TODO: Potentially use the fact that only a very small section of the molecule changes
#at each step in order to optimize segmentation into sources/sinks
def segment_molecule(molecule):
    """
    Segments a molecule into sets of sources and sinks by iterating over all
    source and sink molecular regular expressions this program knows about,
    and creating a data structure that "tags" atoms by their source and sink
    characteristics.

    Returns an object with two properties:

    - sources: List of all sources in the molecule
    - sinks: List of all sinks in the molecule

    Each of these two lists will consist of objects that have the following properties:

    - "subtype" : Type of generic source or sink the object represents, such as "Y-L" or "Y=C-L"
    - "atoms" : An object which has keys for each of the "pieces" of the source/sink,
        and a dict for the actual atom, which has both the atom itself, and its corresponding molecule.
        The last bit exists because Atom objects aren't "real" in that they are generated from scratch
        every time you access a Molecule's atoms.
    """
    result = {"sources" : [], "sinks" : []}
    result["sources"] = label_sources(result["sources"],molecule)
    result["sinks"] = label_sinks(result["sinks"],molecule)
    #do any additional processing here
    return result

def label_sources(sources,molecule):
    """
    Handles the logic of actually looking at the SMILES patterns and producing a list of sources
    based on a molecule.

    Exists to remove clutter from segment_molecule.
    While superficially similar to label_sinks in the basic sense, merging them could get problematic and
    hard to read later.
    """
    for source_type in SOURCES:
        smarts = pybel.Smarts(SOURCES[source_type])
        groups = get_real_indices(smarts.findall(molecule))
        for group in groups:
            source = {"subtype":source_type,"atoms":{}}
            #here is where you wish Python had a switch statement
            if source_type == "Y:":
                #only one atom to label, might as well do it here
                source["atoms"]["Y:"] = {"atom": molecule.atoms[group[0]], "molecule": molecule}
            sources.append(source)
    return sources

def label_sinks(sinks,molecule):
    """
    Handles the logic of actually looking at the SMILES patterns and producing a list of sinks 
    based on a molecule.

    Exists to remove clutter from segment_molecule.
    While superficially similar to label_sources in the basic sense, merging them could get problematic and
    hard to read later.
    """
    for sink_type in SINKS:
        smarts = pybel.Smarts(SINKS[sink_type])
        groups = get_real_indices(smarts.findall(molecule))
        for group in groups:
            sink = {"subtype":sink_type,"atoms":{}}
            if sink_type == "H-L":
                #hydrogens are usually after all the other atoms, but we shouldn't assume this until optimization stage
                for atom_idx in group:
                    if molecule.atoms[atom_idx].atomicnum == 1:
                        sink["atoms"]["H"] = {"atom": molecule.atoms[atom_idx], "molecule": molecule}
                    else:
                        #assume identification was correct, OpenBabel is old.
                        sink["atoms"]["L"] = {"atom": molecule.atoms[atom_idx], "molecule": molecule}
            sinks.append(sink)
    return sinks
