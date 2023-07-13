'''
Carries out the NuL reaction.
'''
from reaction import Reaction
import structure.struct_ops as struct_ops
import structure.scoring as scoring
import structure.properties as properties
import utils
from openbabel import openbabel
from openbabel import pybel

class NuL(Reaction):

    reaction_type = "NuL"
    def __init__(self, sources, sinks):
        Reaction.__init__(self, sources, sinks)
        # For the sink, Y-L, if the two halogens are different, measure the electronegativity of both Halogens. 
        # TODO: Update the code to work on the sink type "ph-S-Cl" 
        self.electrophilic_end = "" # the electrophilic_end of the sink
        self.nucleophilic_end = ""  # The nucleophilic_end of the sink
 
    def cross_check(self):
        if self.cross_check_score != -2.0:
            return self.cross_check_score
        MAGIC_THRESHOLD = 0.3
        self.cross_check_score = MAGIC_THRESHOLD
        # Assign the electrophilic_end and the nucleophilic_end parts of the sink to their crosponding atoms 
        sink = self.sinks[0]
        sink_subtype = sink.subtype
        self.electrophilic_end = "Y" 
        self.nucleophilic_end  = "L"
        mol = sink.molecule
        # If the sink type is Y-L, and the two Halogens are different, measure the electronegativity.
        if mol.OBMol.GetAtom(sink.get_atom("Y")).GetAtomicNum() != mol.OBMol.GetAtom(sink.get_atom("L")).GetAtomicNum():
                electronegativity_Y = OBElements.GetElectroNeg(mol.OBMol.GetAtom(sink.get_atom("Y")).GetAtomicNum())
                electronegativity_L = OBElements.GetElectroNeg(mol.OBMol.GetAtom(sink.get_atom("L")).GetAtomicNum())
                #The Halogen with the largest electronegativity value will be the electrophilic_end, and the other halogen will be the nucleophilic_end
                if electronegativity_Y > electronegativity_L :
                    self.electrophilic_end = "Y" 
                    self.nucleophilic_end = "L" 
                else:
                    self.electrophilic_end = "L"
                    self.nucleophilic_end = "Y" 
        return self.cross_check_score

    def rearrange(self):
		# make "C1-nucleophilic_end" bond and "C2-nucleophilic_end" bond
		# break C=C double bond, break Y-L bond	
		double_bond = self.sources[0] # C=C
		sink = self.sinks[0] # Y-L
		mol = double_bond.molecule
		struct_ops.make_bond(double_bond.get_atom("C1") ,sink.get_atom(self.nucleophilic_end),mol) #make "C1-nucleophilic_end" bond 
		struct_ops.break_bond(sink.get_atom(self.nucleophilic_end) ,sink.get_atom(self.electrophilic_end),mol) # break Y-L bond
		struct_ops.make_bond(sink.get_atom(self.nucleophilic_end), double_bond.get_atom("C2"),mol) # make "C2-nucleophilic_end" bond
		struct_ops.break_bond(double_bond.get_atom("C2") ,double_bond.get_atom("C1"),mol) # break C=C double bond 
