#!/usr/bin/python
#coding=utf-8
"""
Decision-making engine for SiC³. Invoked by sic.py (and by extension sigc.py) to perform mechanistic examinations.
Uses all the segmentation and pka tools in the other modules.
"""

import segmentation.segmentation as segmentation
import pka.pka as pka
import structure.struct_ops as struct_ops
import reaction_types.reaction_factory as reaction_factory #too many of these
import reaction_types.interactions as interactions
import structure.connectivity_table as connectivity_table
import pybel
from reaction_state import ReactionState
import copy
    
def generate_choices(state):
    #first, assign pka and get sources/sinks
    #NOTE: Figure out how to optimize so we don't recalculate these too much, especially pKa
    pka.get_all_pka(state.molecule)
    possible_sites = segmentation.segment_molecule(state.molecule)
    #and now for each source-sink pair, get the interactions
    for source in possible_sites["sources"]:
        for sink in possible_sites["sinks"]:
            interaction_tuple = (source.subtype,sink.subtype)
            possible_interactions = interactions.INTERACTIONS[interaction_tuple] if interaction_tuple in interactions.INTERACTIONS else []
            #listify so that our interface works - if you have multiple sources or multiple sinks, this automagically takes care of it
            #but don't listify if already a list
            interaction_source = source
            interaction_sink = sink
            if type(source) != type([]):
                interaction_source = [source]
            if type(sink) != type([]):
                interaction_sink = [sink]
            for interaction in possible_interactions:
                new_mol = struct_ops.copy_molecule(state.molecule)
                reaction = reaction_factory.produce_reaction(interaction,interaction_source,interaction_sink,mol=new_mol)
                if reaction.cross_check() > 0: #make sure it is actually a possibility
                    new_state = ReactionState(new_mol,parent_state=state,parent_reaction=reaction)
                    #NOTE: A mysterious bug happens wher eif you don't run this line here, suddenly the molecule attached to your sources is not the same as the one on the ReactionState...
                    reaction.rearrange()
                    #DO NOT MOVE REACTION.REARRANGE AWAY FROM HERE
                    state.possibilities.add(new_state)
    #TODO: Make this logging.debug...
    print "%s possibilities" % len(state.possibilities)
    print "cross check of possibilities:"
    for possibility in state.possibilities:
        print possibility.parent_reaction.cross_check()
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
        raise ValueError("No pathh was found between reactants and products.")

def get_mechanism(reactants,products,solvent=False):
    """
    Takes in the reactaants and products as SMILES strings with . separating each molecule in both,
    and returns a list with the ReactionState objects that represent how we got there.
    """
    #TODO: Add check for atom balance, raise ValueError if wrong.
    path_to_product = [] #keeps track of the reaction state that we want to print out at the end
    #path_to_product holds onyl the states that get you to the product, and nothing else in the tree.
    #read in reactants and products
    react_mol = pybel.readstring("smi",reactants)
    prod_mol = pybel.readstring("smi",products)
    react_mol.addh()
    react_mol.connectivity_table = connectivity_table.ConnectivityTable(react_mol)
    prod_mol.addh()
    prod_mol.connectivity_table = connectivity_table.ConnectivityTable(prod_mol)
    #now create a ReactionState out of the reactants - this will be the root
    current_state = ReactionState(react_mol,prod=prod_mol) #product becomes part of the tree
    path_to_product.append(current_state) #since the first state HAS to be the first step in the mechanism
    counter = 0
    print react_mol.write("can")
    while not current_state.matches_product():
        #from current_state, generate choices
        generate_choices(current_state)
        if len(current_state.possibilities) > 0:
            for possibility in current_state.possibilities:
                if not hasattr(possibility,'examined'): #if we didn't look at it and conclude none of its paths get us to product...
                    #check if closer to product
                    if possibility.closer_to_product():
                        current_state = possibility
                        path_to_product.append(possibility)
                        break
            else:
                #none of the possibilities are closer, did not encounter break statement, so go up a level
                current_state = go_up_a_level(current_state,path_to_product)

        else:
            #if there's no further paths and we're still not at product, go up a level too
            current_state = go_up_a_level(current_state,path_to_product)
        counter += 1
        if counter > 20:
            raise ValueError("Could not find a reaction after {} steps.".format(counter))
    #when we hit product, return
    return path_to_product
