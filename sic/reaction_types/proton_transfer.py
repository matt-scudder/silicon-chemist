"""
Carries out the proton transfer reaction.
"""
from .reaction import Reaction
from structure import struct_ops, scoring
from pka import pka

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
            source_sbtype = source.subtype
            # if the sourece type = "Y", get the Pka_BH of "Y"
            # else: source is "Z=C",  get the Pka_BH of "Z"         
            source_atom = "Y" if source_sbtype == "Y" else "Z"
            #get pKa_HA of the H 
            #TODO: Update this code so that it supports more than just Y, maybe use a dict? H-L works always.            
            pKa_HA = pka.get_pka(sink.get_atom("H"),sink.molecule)            
            pKa_BH= pka.get_pka(source.get_atom(source_atom),source.molecule)
            print("pKa_BH =",pKa_BH)
            print("pKa_HA =", pKa_HA)
            if pKa_BH is None or pKa_HA is None:
                #None generally means either "infinite" or "not in our pKa chart".
                #If running debug mode, print(which are None)
#                print("None pKa encountered. pKaBHNu is: {} (atom index {}), pKaHA is: {} (atom index {})".format(pKa_BH,source.get_atom("Y"),pKa_HA,sink.get_atom("H")))
                self.cross_check_score = 0.0
                return self.cross_check_score
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
        source_sbtype = source.subtype
        source_atom = "Y" if source_sbtype == "Y" else "Z"
        #print("source ",source.molecule)
        #TODO: Make this work so it supports more than just Y, see cross_check above.
        #make Y-H bond
        struct_ops.make_bond(source.get_atom(source_atom),sink.get_atom("H"),source.molecule)
        #print("make Y-H bond = ",source.molecule)
        #break H-L bond
        struct_ops.break_bond(sink.get_atom("H"),sink.get_atom("L"),source.molecule)
        if source_sbtype == "Z=C":
            struct_ops.break_bond(source.get_atom("C"),source.get_atom(source_atom),source.molecule)


