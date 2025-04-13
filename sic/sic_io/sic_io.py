"""
Handles all File I/O and related functions.
Primarily used for convenience in sic.py.
"""

import re

from sic.utils import write_mol

SMILES_CHARS = re.escape(r"+[]()=#@/\-") #special characters involved in SMILES, see spec as described in parse_sic_file

def parse_sic_file(sic_input):
    """
    Parses SiC-format file and returns an object with reactants, solvent and products.

    SiC format is defined by the following pattern:

    {reactants}>>{solvent}>>{products}

    With delimiters between each reactant or product being arbitrary, SO LONG AS
    the same delimiter is used between each reactant, product or solvent.
    Characters that CANNOT be used as delimiters are SMILES characters such as =, #,
    [,], etc. (see SMILES format guide at http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html ).

    SiC-format files describe a SINGLE reaction on a SINGLE line. Unlike the original SiC, this program
    does not require a ; character to terminate the line, only EOF or \n; any characters after the last
    SMILES string will be parsed out.

    This function Separates out reactants by taking in all the characters before the first >,
    solvent by taking in any characters found after the first two > in the string,
    and products by taking in all the characters after the last >.
    """

    reaction = sic_input.readlines()
    if len(reaction) > 1:
            raise ValueError("SiC only supports a single reaction at a time.")
    elif len(reaction) < 1:
            raise ValueError("No reactions found in the file provided.")
    reaction_pieces = reaction[0].split(">")
    reactants = reaction_pieces[0]
    products = reaction_pieces[-1].replace(";","").replace("\n","") #some files have more than one \n...
    solvent = False
    if reaction_pieces[2] != "":
            #if you split by ">", then the 2-4th groups will be "" UNLESS there are characters in the middle
            #if there are characters after the second ">", then we parse them in as the solvent
            solvent = reaction_pieces[2]
    react_obj = {"reactants": reactants, "products" : products, "solvent" : solvent}
    for key in react_obj:
            if react_obj[key]:
                    #Since we don't know the arbitrary delimiter, find it using regex and what we know about the SMILES spec
                    react_obj[key] = re.split(r"[^a-zA-Z0-9%s]"%SMILES_CHARS,react_obj[key]) 
    return react_obj


def write_up_mechanism(reaction_state_list,solvent=False):
    """
    Given a list of ReactionState objects, makes a pretty-print representation of the
    mechanism. Similar to SiC's original functionality, but cooler.
        If solvent is provided, inserts it into the first line. Solvent should be added
        as needed into the ReactionStates otherwise.
    Assumes that the list has at least one element, and that the first element is the
    reactants.
    Uses the fact that all ReactionState objects contain a reference to the product.
    """

    reactants_state = reaction_state_list[0]
    reactants = write_mol(reactants_state.molecule)
    product = write_mol(reactants_state.product)
    solv_string = write_mol(solvent) if solvent else ""
    master_string = "%s>>%s>>%s\n"%(reactants,solv_string,product)
    if len(reaction_state_list) < 2:
        master_string += "Trivial mechanism, already at products."
        return master_string
    for i in range(1,len(reaction_state_list)):
        master_string += "Step %s: %s>>%s>>%s\n"%(i,write_mol(reaction_state_list[i-1].molecule),solv_string,write_mol(reaction_state_list[i].molecule))
        master_string += "Reaction type: %s\n" % reaction_state_list[i].parent_reaction.reaction_type
    return master_string

def create_state_smiles(molecules):
    """
    Creates a combined SMILES string from a set of molecular SMILES strings by concatenating them with dots. Nothing fancy.
    """

    combined_smiles = molecules[0]
    for molecule in molecules[1:]:
        combined_smiles += "." + molecule
    return combined_smiles
