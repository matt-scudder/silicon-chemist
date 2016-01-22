#!/usr/bin/python
#coding=utf-8
"""
Decision-making engine for SiC³. Invoked by sic.py (and by extension sigc.py) to perform mechanistic examinations.
Uses all the segmentation and pka tools in the other modules.
"""

from ..segmentation import segmentation
from ..pka import pka
from ..structure import similarity
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
                state.possibility.add(pybel.readstring("smi",state.state.write("smiles")),parent_state=state,parent_reaction=reaction)
    #don't return anything, this just modifies the state and adds in possibilities

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
    current_state = ReactionState(react_mol)
    MASTER_STATE.append(react_mol) #since the first state HAS to be the first step in the mechanism
    while not similarity.is_same_molecule(current_state.state,prod_mol):
        #from current_state, generate choices
        generate_choices(current_state)
        for possibility in current_state.possibilities:
            #rearrange the atoms
            possibility.parent_reaction.rearrange() #move the atoms around
            #check if closer to product
            if similarity.closer_to_product(possibility.state,current_state.state,prod_mol): 
                current_state = possibility
                break
    #when we hit product, return
    return MASTER_STATE
    

    

