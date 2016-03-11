#!/usr/bin/python
#coding=utf-8
"""
Carries out the proton transfer reaction.
"""
from reaction import Reaction
import structure.struct_ops as struct_ops
import pka.pka as pka
import structure.scoring as scoring

class ProtonTransfer(Reaction):
    """
    Reaction object that carries out a proton transfer.
    """
    reaction_type = "proton_transfer"

    def cross_check(self):
        """
        The cross-check score is generated as a simple ΔpKa rule score,
        with -10 or lower being declared extremely unlikely (0.0) and
        10 or higher being declared extremely likely (1.0), with the 0 point
        tilted toward favorable (0.6) and a linear scaling between the points.
        """
        if self.cross_check_score != -2.0: 
            return self.cross_check_score
        else:
            #first get the source and sink in, we only have one so make this readable
            source = self.sources[0]
            sink = self.sinks[0]
            #get pKa_HA of the H, and pKa_BH of the Y
            #TODO: Update this code so that it supports more than just Y, maybe use a dict? H-L works always.
            pKa_BH = pka.get_pka(source.get_atom("Y"),source.molecule)
            pKa_HA = pka.get_pka(sink.get_atom("H"),sink.molecule)
            dpKa = pKa_BH - pKa_HA
            if dpKa < -10: 
                self.cross_check_score = 0.0
            elif dpKa > 10:
                self.cross_check_score = 1.0
            else:
                self.cross_check_score = scoring.score_pka(10,0.6,dpKa)
            return self.cross_check_score

    def rearrange(self):
        """
        Bonds the lone pairs on the basic atom to the H in question,
        and separates the H from the acid.
        """
        source = self.sources[0]
        sink = self.sinks[0]
        #TODO: Make this work so it supports more than just Y, see cross_check above.
        #make Y-H bond
        struct_ops.make_bond(source.get_atom("Y"),sink.get_atom("H"),source.molecule)
        #break H-L bond
        struct_ops.break_bond(sink.get_atom("H"),sink.get_atom("L"),source.molecule)

