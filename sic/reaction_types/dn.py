#!usr/bin/python
#coding=utf-8
"""
Carries out a DN (disassociation of carbon and leaving group)
reaction. Scores based on bonds oon the carbon.
Depends on DUM type for sources because this is essentially
a sink interacting with itself.
"""
from reaction import Reaction
import structure.struct_ops as struct_ops
import structure.properties as properties
import utils
import pka.pka as pka
import structure.scoring as scoring

class DN(Reaction):
    """
    Reaction object that carries out a DN reaction.
    """

    reaction_type = "DN"
    def cross_check(self):
        """
        For a DN, you want to check:
        1. how good the cation is
        2. how good the leaving group is.
        The second check is easily handled by the ΔpKa code from SN2,
        since we need to split off the L to see things and stuff.
        """
        if self.cross_check_score != -2.0:
            return self.cross_check_score
        sink = self.sinks[0]
        FACTORS = {0: 0.0, 1: 0.0, 2: 0.6, 3: 1.0}
        carbon_count = properties.get_carbon_degree(sink)
        final_multiplier = FACTORS[carbon_count]
        if final_multiplier == 0:
            self.cross_check_score = 0.0
            return self.cross_check_score
        #same check for pKa of L as in SN2, but this time there's no ΔpKa, we're just checking L pKa.
        new_mol = struct_ops.copy_molecule(sink.molecule)
        new_sinks = utils.shift_molecule_references(self.sinks,new_mol)
        new_sink = new_sinks[0]
        struct_ops.break_bond(new_sink.get_atom("C"),new_sink.get_atom("L"),new_sink.molecule) #this works because we use indices, not direct atom.OBAtom refs
        #need to break L, C because otherwise C gets a - charge and L gets a + (and an implicit H by SMILES standards...)
        pka.get_all_pka(new_mol)
        pKa_BHL = pka.get_pka(new_sink.get_atom("L"),new_mol)
        if pKa_BHL > 6:
            self.cross_check_score = 0.0
        elif pKa_BHL < -6:
            self.cross_check_score = 1.0
        else:
            self.cross_check_score = scoring.score_pka(6,0.4,dpKa)
        return self.cross_check_score * final_multiplier

    def rearrange(self):
        """
        Breaks the C-L bond.
        """
        sink = self.sinks[0]
        struct_ops.break_bond(sink.get_atom("C"),sink.get_atom("L"),sink.molecule)

