"""
Carries out an SN2 reaction.
"""
from .reaction import Reaction 
from structure import struct_ops, properties, scoring
from pka import pka

class SN2(Reaction):
    """
    Reaction object that carries out an SN2 reaction.
    """
    reaction_type = "SN2"

    def cross_check(self):
        """
        The cross-check score is two-parted.
        One part is the ΔpKa rule: the pKaBH of the leaving group atom
        should be less than 10 units more basic than the incoming nucleophile, or:
        pKaBH(Nu) - pKaBH(L) < -10
        The second part is ranking by softness of the Nu, which is not yet implemented.
        """
        if self.cross_check_score != -2.0:
            return self.cross_check_score
        nucleophile = self.sources[0]
        sink = self.sinks[0]
#        print "Nucleophile atom: {}".format(nucleophile.get_atom("Y"))
#        print "Sink carbon: {}".format(sink.get_atom("C"))
        FACTORS = {0: 1.0, 1: 0.95, 2: 0.6, 3: 0.0} #depends on how many non-carbon and non-L atoms are bound to it
        #first verify that we can even do sn2
        carbon_count = properties.get_carbon_degree(sink) 
        #print "Carbon count: {}".format(carbon_count)
        final_multiplier = FACTORS[carbon_count]
        if final_multiplier == 0:
            self.cross_check_score = 0.0
            return self.cross_check_score #give up before doing slow calculations
        #evaluation of the pKaBH of the Nu is easy - just get it off the atom
        #TODO: Update this code so that it supports more than just Y
        pKa_BHNu = pka.get_pka(nucleophile.get_atom("Y"),nucleophile.molecule)
        #for the L it's harder because the pKa changes based on whether it's bonded to the C
        #we need to break the bond, see how that affects the pKa, and then use that number
        #for now, we copy the molecule, break the bond there, and recalc pKa
        new_mol = struct_ops.copy_molecule(nucleophile.molecule) #the molecule passed in is arbitrary - remember it's the same for source and sink
        struct_ops.break_bond(sink.get_atom("C"),sink.get_atom("L"),new_mol) #this works because we use indices, not direct atom.OBAtom refs
        #need to break L, C because otherwise C gets a - charge and L gets a + (and an implicit H by SMILES standards...)
        pka.get_all_pka(new_mol)
        pKa_BHL = pka.get_pka(sink.get_atom("L"),new_mol)
        #hope the garbage collector kills new_mol here and move on
        if pKa_BHNu is None or pKa_BHL is None:
            #None generally means either "infinite" or "not in our pKa chart".
            #If running debug mode, print which are None
#            print "None pKa encountered. pKaBHNu is: {} (atom index {}), pKaBHL is: {} (atom index {})".format(pKa_BHNu,nucleophile.get_atom("Y"),pKa_BHL,sink.get_atom("L"))
            self.cross_check_score = 0.0
            return self.cross_check_score
        dpKa = pKa_BHNu - pKa_BHL
        if dpKa < -10: 
            self.cross_check_score = 0.0
        elif dpKa > 10:
            self.cross_check_score = 1.0
        else:
            self.cross_check_score = scoring.score_pka(10,0.6,dpKa)
        self.cross_check_score *= final_multiplier
        return self.cross_check_score

    def rearrange(self):
        """
        Breaks the C-L bond and makes the C-Nu bond.
        Stereochem is not taken into account here.
        """
        source = self.sources[0]
        sink = self.sinks[0]
        struct_ops.make_bond(source.get_atom("Y"),sink.get_atom("C"),sink.molecule)
        struct_ops.break_bond(sink.get_atom("C"),sink.get_atom("L"),sink.molecule)
