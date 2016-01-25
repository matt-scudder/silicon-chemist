#!/usr/bin/python
#coding=utf-8
"""
Decision-making engine for SiC³. Invoked by sic.py (and by extension sigc.py) to perform mechanistic examinations.
Uses all the segmentation and pka tools in the other modules.
"""

from ..segmentation import segmentation
from ..pka import pka
from ..structure import struct_ops 
from ..reaction_types import reaction_type_factory #there's going to be too many of these...
from ..reaction_types import interactions
import pybel
from reaction_state import ReactionState
    
def generate_choices(state):
    #first, assign pka and get sources/sinks
    #NOTE: Figure out how to optimize so we don't recalculate these too much, especially pKa
    pka.get_all_pka(state.state)
    possible_sites = segmentation.segment_molecule(state.state)
    #and now for each source-sink pair, get the interactions
    for source in possible_sites["sources"]:
        for sink in possible_sites["sinks"]:
            interaction_tuple = (source["subtype"],sink["subtype"])
            possible_interactions = interactions.INTERACTIONS[interaction_tuple]
            #listify so that our interface works - if you have multiple sources or multiple sinks, this automagically takes care of it
            #but don't listify if already a list
            interaction_source = source
            interaction_sink = sink
            if type(source) != type([]):
                interaction_source = [source]
            if type(sink) != type([]):
                interaction_sink = [sink]
            for interaction in possible_interactions:
                reaction = reaction_type_factory.produce_reaction_type(interaction,interaction_source,interaction_sink)
                if reaction.cross_check() > 0: #make sure it is actually a possibility
                    state.possibility.add(ReactionState(struct_ops.copy_molecule(state.state),parent_state=state,parent_reaction=reaction))
    #don't return anything, this just modifies the state and adds in possibilities

def go_up_a_level(state,master):
    """
    Utility function for get_mechanism.
    Takes a state and goes up a level if there is a level to rise, returning
    the parent of the current state.
    Flags the state as "already examined fully" since this function is called
    when all the possibilities below are known to not get us closer to product, and
    removes it from the master reaction track.
    If there is no level, raises a ValueError with a message indicating the problem.
    """
    if state.parent_state:
        state.examined = True #mark state as examined
        master.pop()
        return state.parent_state
    else:
        raise ValueError("Tried to go up a level but already at root.")

def get_mechanism(reactants,products,solvent=False):
    """
    Takes in the reactaants and products as SMILES strings with . separating each molecule in both,
    and returns a list with the ReactionState objects that represent how we got there.
    """
    MASTER_STATE = [] #keeps track of the reaction state that we want to print out at the end
    #read in reactants and products
    react_mol = pybel.readstring("smi",reactants)
    prod_mol = pybel.readstring("smi",products)
    react_mol.addh()
    prod_mol.addh()
    #now create a ReactionState out of the reactants - this will be the root
    current_state = ReactionState(react_mol,prod=prod_mol) #product becomes part of the tree
    MASTER_STATE.append(react_mol) #since the first state HAS to be the first step in the mechanism
    while not current_state.matches_product():
        #from current_state, generate choices
        generate_choices(current_state)
        if len(current_state.possibilities) > 0:
            for possibility in current_state.possibilities:
                if not possibility.examined: #if we didn't look at it and conclude none of its paths get us to product...
                    #rearrange the atoms
                    possibility.parent_reaction.rearrange() #move the atoms around
                    #check if closer to product
                    if possibility.closer_to_product():
                        current_state = possibility
                        MASTER_STATE.append(possibility)
                        break
            else:
                #none of the possibilities are closer, did not encounter break statement, so go up a level
                current_state = go_up_a_level(current_state,MASTER_STATE)

        else:
            #if there's no further paths and we're still not at product, go up a level too
            current_state = go_up_a_level(current_state,MASTER_STATE)
    #when we hit product, return
    return MASTER_STATE
