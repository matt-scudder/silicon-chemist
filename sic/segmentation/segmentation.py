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
import source
import sink
import pybel

#TODO: Potentially use the fact that only a very small section of the molecule changes
#at each step in order to optimize segmentation into sources/sinks
def segment_molecule(molecule):
    """
    Segments a molecule into sets of sources and sinks by iterating over all
    source and sink molecular regular expressions this program knows about,
    and creating a data structure that "tags" atoms by their source and sink
    characteristics.

    Returns an object with two properties:

    - sources: List of all sources in the molecule, as Source objects
    - sinks: List of all sinks in the molecule, as Sink objects
    """
    result = {"sources" : label_sources(molecule), "sinks" : label_sinks(molecule)}
    #do any additional processing here
    for sink in result["sinks"]:
        #check for C-Ls that might fall apart
        if sink.subtype == "C-L":
            #add a dummy source to sources, just one is enough to pair with all the C-Ls
            result["sources"].append(source.Source("DUM",[],molecule))
            break
    return result

def label_sources(molecule):
    """
    Handles the logic of actually looking at the SMILES patterns and producing a list of sources
    based on a molecule.

    Exists to remove clutter from segment_molecule.
    While superficially similar to label_sinks in the basic sense, merging them could get problematic and
    hard to read later.
    """
    sources = []
    mol_atoms = molecule.atoms #so that we don't do a list processing every time, given that molecule.atoms would regenerate itself
    for source_type in SOURCES:
        smarts = pybel.Smarts(SOURCES[source_type])
        groups = smarts.findall(molecule)
        for group in groups:
            sources.append(source.Source(source_type,group,molecule))
    return sources

def label_sinks(molecule):
    """
    Handles the logic of actually looking at the SMILES patterns and producing a list of sinks 
    based on a molecule.

    Exists to remove clutter from segment_molecule.
    While superficially similar to label_sources in the basic sense, merging them could get problematic and
    hard to read later.
    """
    sinks = []
    mol_atoms = molecule.atoms
    for sink_type in SINKS:
        smarts = pybel.Smarts(SINKS[sink_type])
        groups = smarts.findall(molecule)
        for group in groups:
            sinks.append(sink.Sink(sink_type,group,molecule))
    return sinks
