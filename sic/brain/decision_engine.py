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
import structure.properties as properties
import utils
from openbabel import pybel
from reaction_state import ReactionState
import copy
import time
    
def generate_choices(state):
    #first, assign pka and get sources/sinks
    #NOTE: Figure out how to optimize so we don't recalculate these too much, especially pKa
    pka.get_all_pka(state.molecule)
#    print "Current state: {}".format(state.molecule.write("can"))
#    print "Initial bond distance: {}".format(properties.get_bond_distance(state.molecule,state.product,state.mapping)) 
    possible_sites = segmentation.segment_molecule(state.molecule)
    #and now for each source-sink pair, get the interactions
    for source in possible_sites["sources"]:
        for sink in possible_sites["sinks"]:
            interaction_tuple = (source.subtype,sink.subtype)
#            print "Interaction: {}".format(interaction_tuple)
            possible_interactions = interactions.INTERACTIONS[interaction_tuple] if interaction_tuple in interactions.INTERACTIONS else []
            #listify so that our interface works - if you have multiple sources or multiple sinks, this automagically takes care of it
            #but don't listify if already a list
#            print "Possible interactions: {}".format(possible_interactions)
            interaction_source = source
            interaction_sink = sink
            if type(source) != type([]):
                interaction_source = [source]
            if type(sink) != type([]):
                interaction_sink = [sink]
            for interaction in possible_interactions:  
                new_mol = struct_ops.copy_molecule(state.molecule)
                reaction = reaction_factory.produce_reaction(interaction,interaction_source,interaction_sink,mol=new_mol)
#                print "Made reaction of type {}, cross check is {}".format(interaction,reaction.cross_check())
                if reaction.cross_check() > 0: #make sure it is actually a possibility
                #TODO: Add a method to the Reaction class that tells whether there are two copies of the source or the sink to check
                    new_state = ReactionState(new_mol,parent_state=state,parent_reaction=reaction)
                    #NOTE: A mysterious bug happens where if you don't run this line here, suddenly the molecule attached to your sources is not the same as the one on the ReactionState...
                    reaction.rearrange()
                    #DO NOT MOVE REACTION.REARRANGE AWAY FROM HERE
                    state.possibilities.add(new_state)
                    two_products_interactions = ["AE","ADE3","E2"]#Tif the interaction is "AE","ADE3", or "E2", there should be two correct products 
                    if interaction in two_products_interactions:
                        new_mol = struct_ops.copy_molecule(state.molecule)
                        reaction = reaction_factory.produce_reaction(interaction,interaction_source,interaction_sink,mol=new_mol, second_product = True)
                        print "Made reaction2 of type {}, cross check is {}".format(interaction,reaction.cross_check())
                        if reaction.cross_check() > 0:
                            new_state = ReactionState(new_mol,parent_state=state,parent_reaction=reaction)
                            reaction.rearrange()
                            state.possibilities.add(new_state)
    #TODO: Make this logging.debug...
#    print "%s possibilities" % len(state.possibilities)
#    print "Current state is: %s" % state.molecule.write("can")
#    print "cross check of possibilities:"
#    for possibility in state.possibilities:
#        print possibility.parent_reaction.reaction_type
#        print possibility.parent_reaction.cross_check()
#        print possibility.molecule.write("can")
#        print "Bond distance: {}, Closer to product: {}".format(properties.get_bond_distance(possibility.molecule,possibility.product,possibility.mapping),possibility.closer_to_product())
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

def get_mechanism(reactants,products,solvent=False,max_counter=15):
    """
    Takes in the reactaants and products as SMILES strings with . separating each molecule in both,
    and returns a list with the ReactionState objects that represent how we got there.
    """
    start_time = time.time()
    #TODO: Add check for atom balance, raise ValueError if wrong.
    path_to_product = [] #keeps track of the reaction state that we want to print out at the end
    #path_to_product holds onyl the states that get you to the product, and nothing else in the tree.
    #read in reactants and products
    react_mol = pybel.readstring("smi",reactants)
    prod_mol = pybel.readstring("smi",products)
    #now create a ReactionState out of the reactants - this will be the root
    before_java = time.time()
    current_state = ReactionState(react_mol,prod=prod_mol) #product becomes part of the tree
    after_java = time.time()
    #doing the ReactionState initializer changes the internal reactant and product, so put them back in
    #TODO: Make this cleaner
    react_mol = current_state.molecule
    prod_mol = current_state.product
    react_mol.removeh() #this looks stupid, but sometimes hydrogens are added explicitly, counteracting our assumption that all backbone atoms come before all H atoms
    #since breaking this assumption makes bond distance stop working, this seemingly-stupid function call is VITAL and should NOT BE REMOVED
    react_mol.addh()
    prod_mol.removeh()
    prod_mol.addh()
    react_mol.connectivity_table = connectivity_table.ConnectivityTable(react_mol)
    prod_mol.connectivity_table = connectivity_table.ConnectivityTable(prod_mol)
    print react_mol.connectivity_table.closer_to_product_table
    path_to_product.append(current_state) #since the first state HAS to be the first step in the mechanism
    counter = 0
    print react_mol.write("can")
    print prod_mol.write("can")
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
        if counter >= max_counter:
            raise ValueError("Could not find a reaction after {} steps.".format(counter))
    #when we hit product, return
    final_time = time.time()
    print "Total time (with java): %s" % (final_time - start_time)
    print "Total time (without java): %s" % (final_time - start_time - (after_java - before_java))
    return path_to_product
