#!/usr/bin/python
#coding=utf-8
"""
Base class that all reaction types extend.
Defines basic methods that all reaction types
should be able to execute, as well as data members.
"""

class ReactionType():
    """
    A Reaction is defined as an interaction between a list of sources
    and a list of sinks. While usually only one source and one sink
    will be involved, this is left open for strange reaction types.

    A Reaction has a cross check, and a reversible rearrangement.

    A cross check is defined as a set of checks
    on the interaction between the source and the sink which produce
    a score between 0 and 1, usually involving the ΔpKa rule.

    A rearrangement is defined as the shift in bonds between the source
    and the sink, resulting in a new molecular state. It is referred to as
    "reversible", because the object should be able to reverse the bond changes
    made in the case that this reaction path does not get us closer to product.

    Because of the fact that we need to keep discrete molecular state in order
    to list the final reaction pathway, as well as walk back and forth between it,
    a ReactionType should take as input a #copy# of the current molecular state,
    and to leave the state as close to it originally was if calling.
    """
    def __init__(self, sources, sinks):
        """
        Attaches sources and sinks to this object, and sets the initial cross check
        score to a sentinel value, so that cross_check() needs only to be run once.
        """
        self.sources = sources
        self.sinks = sinks
        self.cross_check_score = -2.0 

    def cross_check():
        """
        Produces a score on a [0,1] scale based on how likely the source-sink
        interaction attached to this reaction is. Exceedingly unlikely (e.g. ΔpKa < -10)
        reactions should be given a score of 0, and exceedingly likely (e.g. ΔpKa > 10)
        reactions should be given a score of 1, though the latter is subject to change.
        """
        raise NotImplementedError("This is an abstract class and should not be used.")

    def rearrange():
        """
        Rearranges the atoms in the source and sink, in effect carrying out the reaction
        by making the connectivity changes necessary.
        """
        raise NotImplementedError("This is an abstract class and should not be used.")

    def disarrange():
        """
        Reverses anything done by rearrange. It might be necessary to put in a boolean
        that prevents this method from being called when rearrange has not been,
        but for now the assumption is that the convention will prevent callers
        from doing this.
        """
        raise NotImplementedError("This is an abstract class and should not be used.")
