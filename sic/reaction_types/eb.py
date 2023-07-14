"""
Carries out a beta elimination (second step of E1cb).
"""

from .reaction import Reaction 
from sic.structure import struct_ops

class EB(Reaction):
    """
    Reaction object that carries out a beta elimination reaction.
    """

    reaction_type = "EB"

    def cross_check(self):
        """
        An EB (for now) only happens between a [C-] and a C-L.
        Given that the C- can only be produced by a preceding proton_transfer, which itself has a good cross-check
        for ewg and similar, we don't need to do a cross-check here and can just say it will happen,
        like with AN, especially since there is no adjacent C-H to check.
        """

        if self.cross_check_score != -2.0:
            return self.cross_check_score
        self.cross_check_score = 1.0
        return self.cross_check_score

    def rearrange(self):
        """
        Breaks C-L and makes C-C double bond.
        """

        base = self.sources[0] #C-
        CL = self.sinks[0]
        struct_ops.make_bond(base.get_atom("C-"),CL.get_atom("C"),CL.molecule)
        struct_ops.break_bond(CL.get_atom("C"),CL.get_atom("L"),CL.molecule)
