"""
Carries out the "ADN" reaction.
"""
from .reaction import Reaction
from pka import pka
from structure import scoring, struct_ops

class ADN(Reaction):

    reaction_type = "ADN"
    def __init__(self, sources, sinks):
        Reaction.__init__(self, sources, sinks)
        # Since we have two types of sinks,"Z=C" and "C=C",
        # "self.F_atom" will be the first atom  of the double bond; eaither "Z" or "C1"
        # "self.S_atom" will be the second atom  of the double bond; eaither "C" or "C2" 
        # "self.F_atom" and "self.S_atom"  will be used in making and breaking bonds methods
        self.F_atom  = ""
        self.S_atom = ""
 
    def cross_check(self):
        """
        for ADN, Check the ΔpKa rule: :
        pKaBH of the formed carbanin atom should be less than 10 units more basic than the incoming nucleophile, or:
        pKaBH(Nu) - pKaBH(C-) < -10           
        """
        if self.cross_check_score != -2.0 :
           return self.cross_check_score
        nucleophile = self.sources[0]
        sink = self.sinks[0]
        #evaluation of the pKaBH of the Nu is easy - just get it off the atom
        source_sbtype = nucleophile.subtype
        sink_sbtype = sink.subtype
        self.F_atom = "Z" if sink.subtype == "Z=C" else "C1"
        self.S_atom = "C" if sink.subtype == "Z=C" else  "C2"
        pKa_BHNu = pka.get_pka(nucleophile.get_atom(source_sbtype),nucleophile.molecule)
        #for the ewg-C1 it's harder because the pKa changes based on whether it's bonded to the C2
        #we need to break the bond C1=C2 double bond
        #for now, we copy the molecule, break the bond there
        # use a smart pattern to find "ewg-C-" part of the molecule, and recalc pKa
        # for the Z=C, break the double bond, make a C-Nu bond 
        new_mol = struct_ops.copy_molecule(nucleophile.molecule)
        struct_ops.break_bond(sink.get_atom(self.S_atom),sink.get_atom(self.F_atom),new_mol)
        pka.get_all_pka(new_mol)
        pKa_BH = pka.get_pka(sink.get_atom(self.F_atom),new_mol) # get
        if pKa_BHNu is None or pKa_BH is None:
            self.cross_check_score = 0.0
            return self.cross_check_score
        dpKa = pKa_BHNu - pKa_BH
        if dpKa < -10: 
            self.cross_check_score = 0.0
        elif dpKa > 10:
            self.cross_check_score = 1.0
        else:
            self.cross_check_score = scoring.score_pka(10,0.6,dpKa)
        return self.cross_check_score
       
    def rearrange(self):
        # break the double bond , make C bonds to the lone pairs on the source atom
        nucleophile = self.sources[0] # Y, C-
        sink = self.sinks[0]     # ewg-C1=C2 or "Z=C"
        source_sbtype = nucleophile.subtype
        struct_ops.make_bond(nucleophile.get_atom(source_sbtype),sink.get_atom(self.S_atom),sink.molecule)
        struct_ops.break_bond(sink.get_atom(self.S_atom),sink.get_atom(self.F_atom),sink.molecule)