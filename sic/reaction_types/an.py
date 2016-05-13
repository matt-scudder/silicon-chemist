#!/usr/bin/python
#coding=utf-8
"""
Carries out an AN (association to positively-charged carbon center)
reaction. Currently is just a 1.0 score, since
rearrangements are not implemented.
"""
from reaction import Reaction
import structure.struct_ops as struct_ops

class AN(Reaction):
    """
    Reaction object that carries out an AN reaction.
    """
    reaction_type = "AN"

    def cross_check(self):
        """
        Returns 1.0 until rearrangements are implemented.
        """
        if self.cross_check_score != -2.0:
            return self.cross_check_score
        self.cross_check_score = 1.0
        return self.cross_check_score

    def rearrange(self):
        """
        Attaches the nucleophile to the C+.
        """
        source = self.sources[0]
        sink = self.sinks[0]
        #TODO: Update for nucleophiles beyond Y
        struct_ops.make_bond(source.get_atom("Y"),sink.get_atom("C+"),sink.molecule)
        #no bond-breaking involved in an AN.
