#!/usr/bin/python
#coding=utf-8
"""
Carries out an E2 reaction.
"""
from reaction import Reaction 
import structure.struct_ops as struct_ops
import structure.properties as properties

class E2(Reaction):
    """
    Reaction object that carries out an E2 reaction.
    """
    reaction_type = "E2"

    def cross_check(self):
        """
        An E2 (for now) only happens between a lone pair source (Y) and a C-L.
        This cross-check only consists of checking whether C-H is next to C-L.
        While we could check for the whole substitution vs. elimination pathway here,
        remember that we know what the product is, so this path will be discarded if it
        turns out we haven't actually made bonds present in the product.
        That said, it is possible to implement it by checking the degree of the C-L
        carbon, and associating a score multiplier.
        Note that this does not check sterics on the base.
        """
        if self.cross_check_score != -2.0:
            return self.cross_check_score
        base = self.sources[0]
        #NOTE: Maybe put some pKa_BH checks on the base here??
        CL = self.sinks[0]
        if properties.get_adjacent_ch(CL):
            self.cross_check_score = 1.0
        else:
            self.cross_check_score = 0.0
        return self.cross_check_score

    def rearrange(self):
        """
        Breaks C-L and adjacent C-H bonds, makes C-C double bond.
        """
        base = self.sources[0]
        CL = self.sinks[0]
        mol = CL.molecule #because this gets typed a lot
        C_to_change = properties.get_adjacent_ch(CL)
        #get the H, doesn't matter which because we don't care about stereochem yet
        H_to_change = properties.get_H_bonds(C_to_change,CL.mol)[0]
        struct_ops.make_bond(base.get_atom("Y"),H_to_change,mol)
        struct_ops.break_bond(H_to_change,C_to_change,mol)
        struct_ops.break_bond(CL.get_atom("C"),CL.get_atom("L"),mol)
        #remember to make bond from a SOURCE of electrons (i.e. the C that just had the C-H broken)
        #to a SINK of electrons (i.e. the C-L that just had the L leave)
        struct_ops.make_bond(C_to_change,CL.get_atom("C"),mol)
