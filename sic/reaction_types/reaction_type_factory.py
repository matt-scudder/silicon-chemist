#!/usr/bin/python
#coding=utf-8
"""
Factory pattern function that just gives out the right ReactionType given an input string.
This is done because extending the system to take into account all the reaction types as we add them
could get extremely annoying.

NOTE: The reaction types given MUST be the same as the ones in interactions.py!
"""

from proton_transfer import ProtonTransfer

def produce_reaction_type(r_type,sources,sinks):
    """
    Returns the correct reaction type given sources and sinks.
    """
    if r_type == "proton_transfer":
        return ProtonTransfer(sources,sinks)
