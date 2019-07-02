'''
Carries out the Ae, Electrophile Addition to a Multiple Bond, reaction.
'''
from reaction import Reaction
import structure.struct_ops as struct_ops
import structure.scoring as scoring
import structure.properties as properties
import utils
from openbabel import OBElementTable

class AE(Reaction):

    reaction_type = "AE"
    def __init__(self, sources, sinks):
    	Reaction.__init__(self, sources, sinks)
    	# For the source "C=C", mark should be the carbon with the heighest carbon degree on one side ofthe double bond, and "ant_mark" is the carbon that has less smaller carbon degree
        self.mark  = ""
        self.anti_mark = ""
        #Since there are two types of sinks, H-L and Y-L, we will define the electrophilic_end of the sink to be eaither H or Y depending on the sink type
        # The nucleophilic_end of the sink will be L
        # In case the sink type is Y-L, if the two halogens are different, measure the electronegativity of both Halogens. 
        self.electrophilic_end = ""
        self.nucleophilic_end = ""

    def cross_check(self):
        '''
        for AE, check:
            The Mark. rule for each carbon on the edges of the double bond to check"The stability of the Carboncation". 
        	Then compare the two stabilities and determine which carbon will be the "Carboncation" if the crosscheck is not 0.  
        '''
        if self.cross_check_score != -2.0 :
           return self.cross_check_score
        source = self.sources[0]
        FACTORS = {0: 0.0, 1: 0.0, 2: 0.25, 3: 1.0}
        # Getting the carbond degree of the first carbon in the double bond "C1=C2"
        self.carbon_count1 = properties.get_carbon_degree(source,"C1")
        final_multiplier1 = FACTORS[self.carbon_count1]
        # Getting the carbond degree of the first carbon in the double bond
        self.carbon_count2 = properties.get_carbon_degree(source,"C2")
        final_multiplier2 = FACTORS[self.carbon_count2]
        if final_multiplier1 == 0 and final_multiplier2 == 0 :
           self.cross_check_score = 0.0
           return self.cross_check_score
        elif final_multiplier1 > final_multiplier2 :
            self.cross_check_score = final_multiplier1
            self.mark = "C1"
            self.anti_mark = "C2"
            return self.cross_check_score
        else:
            self.cross_check_score = final_multiplier2
            self.mark = "C2"
            self.anti_mark = "C1"           
        # Assign the electrophilic_end and the nucleophilic_end parts of the sink to their crosponding atoms 
        sink = self.sinks[0]
        sink_subtype = sink.subtype
        self.electrophilic_end = "Y" if sink_subtype == "Y-L" else "H"
        self.nucleophilic_end  = "L"
        mol = sink.molecule
        # If the sink type is Y-L, and the two Halogens are different, measure the electronegativity.
        if sink_subtype == "Y-L" and mol.OBMol.GetAtom(sink.get_atom("Y")).GetAtomicNum() != mol.OBMol.GetAtom(sink.get_atom("L")).GetAtomicNum():
                periodic_table = OBElementTable()
                electronegativity_Y = periodic_table.GetElectroNeg(mol.OBMol.GetAtom(sink.get_atom("Y")).GetAtomicNum())
                electronegativity_L = periodic_table.GetElectroNeg(mol.OBMol.GetAtom(sink.get_atom("L")).GetAtomicNum())
                #The Halogen with the largest electronegativity value will be the electrophilic_end, and the other halogen will be the nucleophilic_end
                if electronegativity_Y > electronegativity_L :
                    self.electrophilic_end = "Y" 
                    self.nucleophilic_end = "L" 
                else:
                    self.electrophilic_end = "L"
                    self.nucleophilic_end = "Y" 
        return self.cross_check_score

    def rearrange(self):
    	# make C-H or C-Y bond ,and adjacent C-L bond
    	# break C=C double bond and the H-L bond or the Y-L
    	double_bond = self.sources[0] # C=C
    	sink = self.sinks[0] # H-L for now, update it for more sinks in the future
    	mol = double_bond.molecule
    	struct_ops.break_bond(double_bond.get_atom(self.mark),double_bond.get_atom(self.anti_mark),mol) # break C=C double bond	
    	struct_ops.break_bond(sink.get_atom(self.nucleophilic_end),sink.get_atom(self.electrophilic_end),mol) # 
    	struct_ops.make_bond(double_bond.get_atom(self.anti_mark) ,sink.get_atom(self.nucleophilic_end),mol) 
