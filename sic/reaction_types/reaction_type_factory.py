#!/usr/bin/python
#coding=utf-8
"""
Factory pattern function that just gives out the right ReactionType given an input string.
This is done because extending the system to take into account all the reaction types as we add them
could get extremely annoying.

NOTE: The reaction types given MUST be the same as the ones in interactions.py!
"""

from proton_transfer import ProtonTransfer
from structure.struct_ops import copy_molecule
from pybel import Molecule
import utils


def produce_reaction_type(r_type,sources,sinks,mol=False):
    """
    Returns the correct reaction type given sources and sinks.
    If a molecule is provided, copies the information in sources and sinks
    to be attached to the new molecule.
    """
    orig_mol =  sources[0]["atoms"][sources[0]["atoms"].keys()[0]]["molecule"]
    if mol:
        new_sources = utils.shift_molecule_references(sources,mol)
        new_sinks = utils.shift_molecule_references(sinks,mol)
    else:
        new_sources = sources
        new_sinks = sinks
    if r_type == "proton_transfer":
        return ProtonTransfer(new_sources,new_sinks)
