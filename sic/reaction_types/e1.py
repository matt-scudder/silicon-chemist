"""
Carries out the second step of an E1 elimination (i.e. the elimination part).
"""

from .reaction import Reaction 
from sic.structure import struct_ops, properties

class E1(Reaction):
    """
    Reaction object that carries out the second step of an E1 elimination.
    """
    
    reaction_type = "E1"

    def cross_check(self):
        """
        An E1 (for now) only happens between a [C+] and a Y.
        Given that the [C+] can only be produced by a DN, which itself has a heavy cross-check, we
        only need to check whether there's a C-H to remove here.
        """

        if self.cross_check_score != -2.0:
            return self.cross_check_score
        C = self.sinks[0] #C+
        if properties.get_adjacent_ch(C,carbon_label="C+"):
            self.cross_check_score = 1.0
        else:
            self.cross_check_score = 0.0
        return self.cross_check_score

    def rearrange(self):
        """
        Breaks C-H next to C+, makes Y-H and C+ - C-H bond.
        """

        base = self.sources[0] #Y
        cation = self.sinks[0] #C+
        mol = cation.molecule
        C_to_change = properties.get_adjacent_ch(cation,carbon_label="C+") #should exist since if cross_check_score == 0 then we won't call this method
        H_to_change = properties.get_H_bonds(C_to_change,mol)[0]
        struct_ops.make_bond(base.get_atom("Y"),H_to_change,mol)
        struct_ops.make_bond(C_to_change,cation.get_atom("C+"),mol)
        struct_ops.break_bond(H_to_change,C_to_change,mol)
