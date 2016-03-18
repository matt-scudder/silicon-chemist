#!/usr/bin/python
#coding=utf-8
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
import json
import sortedcontainers
import structure.similarity as similarity
import structure.properties as properties

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
        if not type(self).product and prod:
            type(self).product = prod
        if not type(self).mapping:
            #initialize mapping, do NOT redo mapping, ever!
            type(self).mapping = properties.get_mapping(self.molecule,type(self).product)

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
	json_dict["state"] = molecule.write("can") #this state need not be preserved up top
	json_dict["possibilities"] = map(lambda child: child.to_json_dict(), self.possibilities)
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
        If called on the root node, this method raises a ValueError.
        """
        if not self.parent_state:
            raise ValueError('Root node cannot be checked for "closer to product" since it is the base.')
        return properties.get_bond_distance(self.molecule,self.product,self.mapping) < properties.get_bond_distance(self.parent_state.molecule,self.product,self.mapping)
