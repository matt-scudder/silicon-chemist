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
    - The Molecule object that represents this reaction state (called "state").

At any point where molecules need to be rearranged by a mechanism, a new Molecule should be created by doing mol.write("smiles") and using readstring
on that SMILES string. This new Molecule should then be rearranged before creating a new ReactionState.
"""
import json

class ReactionState():
    def __init__(self,molecule,parent_state=None,parent_reaction=None):
        self.state = molecule
        self.parent_state = None #doesn't matter if None gets assigned
        self.parent_reaction = None
        self.possibilities = []
    
    def to_json_dict(self):
        """
	Writes out a dictionary that will be equivalent to a JSON representation of this object,
	and will be regenerated in the browser by traversing the tree in order to make the children
	have parent references.
	"""
	json_dict = {}
	json_dict["parent_reaction"] = self.parent_reaction.to_json_dict() 
	json_dict["state"] = molecule.write("smiles") #this state need not be preserved up top
	json_dict["possibilities"] = map(lambda child: child.to_json_dict(), self.possibilities)
	return json_dict

