#!/usr/bin/python
#coding=utf-8
"""
Carries out the proton transfer reaction.
"""
from reaction_type import ReactionType
from ..structure import struct_ops
from ..pka import pka

class ProtonTransfer(ReactionType):
    """
    ReactionType object that carries out a proton transfer.
    """

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
	    pKa_BH = pka.get_pka(source["atoms"]["Y:"]["atom"],source["atoms"]["Y:"]["molecule"])
            pKa_HA = pka.get_pka(sink["atoms"]["H"]["atom"],sink["atoms"]["H"]["molecule"])
            dpKa = pKa_BH - pKa_HA
            if dpKa < -10: 
                self.cross_check_score = 0.0
                return self.cross_check_score
            elif dpKa <= 0:
                self.cross_check_score = 0.6 - (0.6*abs(dpKa/ 10.0))
                return self.cross_check_score
            elif dpKa > 10:
                self.cross_check_score = 1.0
                return self.cross_check_score
            else:
                self.cross_check_score = 0.6 + (0.4*abs(dpKa/ 10.0))
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
        struct_ops.make_bond(source["atoms"]["Y:"],sink["atoms"]["H"])
        #break H-L bond
        struct_ops.break_bond(sink["atoms"]["L"],sink["atoms"]["H"])

