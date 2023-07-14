'''
Carries out the Ade3 reaction.
'''
from .reaction import Reaction
from structure import struct_ops, properties

class ADE3(Reaction):

    reaction_type = "ADE3"
    def __init__(self, sources, sinks,second_product = False ):
        Reaction.__init__(self, sources, sinks,second_product)
        # "mark" should be the carbon with the heighest carbon degree on one side ofthe double bond, and "ant_mark" is the carbon that has less smaller carbon degree
        self.mark  = ""
        self.anti_mark = ""

    def cross_check(self):
        '''
        for ADE3, we will assign a threshold to "cross_check_score", so if the "AE" and "ADN" reactions don't pass the cross_check, we prefer "ADE3"
        We still have to check Mark rule for the rearangment function, to decide which of the C1=C2 is mark and anti_mark
        '''
        print("self.second_product =",self.second_product)
        print(" ")
        if self.cross_check_score != -2.0:
            return self.cross_check_score
        source = self.sources[0]
        source_sbtype = source.subtype # type = "C1=C2" or "Z=C"
        if source_sbtype == "C=C" : # if the source is "C1=C2", determin which is mark and anti-mark carbons. 
            FACTORS = {0: 0.0, 1: 0.0, 2: 0.25, 3: 1.0}
            # Getting the carbond degree of the first carbon in the double bond "C1=C2"
            self.carbon_count1 = properties.get_carbon_degree(source,"C1")
            final_multiplier1 = FACTORS[self.carbon_count1]
            # Getting the carbond degree of the second carbon in the double bond
            self.carbon_count2 = properties.get_carbon_degree(source,"C2")
            final_multiplier2 = FACTORS[self.carbon_count2]
            if final_multiplier1 == 0 and final_multiplier2 == 0 :
                self.cross_check_score = 0.                                    
            elif (final_multiplier1 > final_multiplier2) or ((final_multiplier1 == final_multiplier2) and (self.second_product == True)): 
                self.cross_check_score = final_multiplier1
                print("self.cross_check_score2 ADE3 =", self.cross_check_score)
                self.mark = "C1"
                self.anti_mark = "C2"                     
            else:
                self.cross_check_score = final_multiplier2
                print("self.cross_check_score3 ADE3 =", self.cross_check_score)
                self.mark = "C2"
                self.anti_mark = "C1"                    
        MAGIC_THRESHOLD = 0.3
        self.cross_check_score = MAGIC_THRESHOLD
        return self.cross_check_score

    def rearrange(self):
        # make C-H or Z-H bond and adjacent C-L bond
        # break C=C or Z=C double bond and H-L bond
        double_bond = self.sources[0] 
        H_L = self.sinks[0] # H-L for now, update it for more sinks in the future
        mol = double_bond.molecule
        source_sbtype = double_bond.subtype # C=C or Z=C
        self.mark = "C" if source_sbtype == "Z=C" else self.mark
        self.anti_mark = "Z" if source_sbtype == "Z=C" else self.anti_mark
        struct_ops.make_bond(double_bond.get_atom(self.mark) ,H_L.get_atom("L") ,mol) # C-L 
        struct_ops.break_bond(double_bond.get_atom(self.anti_mark) , double_bond.get_atom(self.mark),mol) # break C=C or Z=C double bond
        struct_ops.make_bond(H_L.get_atom("H") ,double_bond.get_atom(self.anti_mark),mol) # C-H bond
        struct_ops.break_bond( H_L.get_atom("L"), H_L.get_atom("H"),mol) # break H-L double bond
        