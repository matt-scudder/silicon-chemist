"""
ReactionState forms the tree-like structure that holds all of the states involved in the reaction.
The entire tree needs to be kept in order to later on tell students whether their choice was the same as they did,
and point out whether it's a reasonable step without having to go through the entire mechanistic shenanigans.

Each ReactionState keeps track of:
    - The state that got it there (if any; root has none)
    - The ReactionType that got it there (if any; root has none). This contains the cross_check score, which allows you to give the user messages.
    - The possible paths from this ReactionState, which are themselves ReactionStates.
    - The Molecule object that represents this reaction state.

All ReactionStates share a "product" object, which ensures the tree structure and allows for utility functions
related to how close a ReactionState is to product. Modify self.product at your own risk.

At any point where molecules need to be rearranged by a mechanism, a new Molecule should be created by doing mol.write("can") and using readstring
on that SMILES string. This new Molecule should then be rearranged before creating a new ReactionState.

This class also contains utility functions to check whether a reaction state is closer to product than its parent,
and whether it matches the product exactly.
"""

from openbabel.pybel import readstring
import sortedcontainers

from sic.structure import similarity, properties

class ReactionState(object):
    product = None
    mapping = None #one mapping from reactants to products for all states, because it is a property of each atom...
    def __init__(self,molecule,parent_state=None,parent_reaction=None,prod=None):
        self.molecule = molecule
        self.parent_state = parent_state #doesn't matter if None gets assigned
        self.parent_reaction = parent_reaction
        #sortedcontainers sorts from lowest to highest, so we sort it "backwards".
        #since cross_check() is on [0,1], subtract score from 1 and take abs.
        #in this scheme, 1 -> 0 and 0 -> 1, making it go in the right order.
        self.possibilities = sortedcontainers.SortedListWithKey(key=lambda x: 1.0 - x.parent_reaction.cross_check())
        if prod:
            #These are sort of static but not really. We want to be able to replace these as we go,
            #because otherwise when this runs as a webserver we have trouble.
            type(self).product = prod
            #initialize mapping, do NOT redo mapping, ever!
            #overwrite product and current_state when given mapping
            #no separate if clause for self.product since mapping and product should be defined at the same time
            mapping,really_canonical_reactants,really_canonical_products = properties.get_mapping(self.molecule,type(self).product)
            type(self).mapping = mapping
            type(self).product = readstring("smi",really_canonical_products)
            #overwrite current state as well
            self.molecule = readstring("smi",really_canonical_reactants)

    def __repr__(self):
        """
        Returns a representation of this object for easy debugging.
        """
        return "ReactionState<ParentReaction:{},Possibilities:{}>".format(self.parent_reaction,self.possibilities)
    
    def to_json_dict(self):
        """
        Writes out a dictionary that will be equivalent to a JSON representation of this object,
        and will be regenerated in the browser by traversing the tree in order to make the children
        have parent references.
        """
        json_dict = {}
        json_dict["parent_reaction"] = self.parent_reaction.to_json_dict() 
        json_dict["state"] = self.molecule.write("can") #this state need not be preserved up top
        json_dict["possibilities"] = [child.to_json_dict() for child in self.possibilities]
        return json_dict

    def matches_product(self):
        """
        Checks whether this reaction state is equal to the product.
        """
        return similarity.is_same_molecule(self.molecule,self.product)

    def closer_to_product(self):
        """
        Examines how close the current ReactionState is to product, and returns True if it is
        closer than its parent.
        If the ReactionState matches product, it will return True, overriding the mapping.
        This is done to avoid the cases where the mapping is wrong.
        If called on the root node, this method raises a ValueError.
        """
        if not self.parent_state:
            raise ValueError('Root node cannot be checked for "closer to product" since it is the base.')
        if self.matches_product():
            return True #circumvent mapper issues
        #This used to be one line, but there's strange failure cases, so let's track them
        current_distance = properties.get_bond_distance(self.molecule,self.product,self.mapping)
        previous_distance = properties.get_bond_distance(self.parent_state.molecule,self.product,self.mapping)
        if current_distance == 0 and not self.matches_product(): #mapper or frozenset issue, investigate this later
            return False
        return current_distance < previous_distance
