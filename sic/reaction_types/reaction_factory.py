"""
Factory pattern function that just gives out the right Reaction given an input string.
This is done because extending the system to take into account all the reaction types as we add them
could get extremely annoying.

NOTE: The reaction types given MUST be the same as the ones in interactions.py!
"""

from .reaction_ade3 import ADE3
from .reaction_adn import ADN
from .reaction_ae import AE
from .reaction_an import AN
from .reaction_dn import DN
from .reaction_e1 import E1
from .reaction_e2 import E2
from .reaction_eb import EB
from .reaction_NuL import NuL
from .reaction_proton_transfer import ProtonTransfer
from .reaction_sn2 import SN2
from sic import utils

def produce_reaction(r_type,sources,sinks,mol=False,second_product = False):
    """
    Returns the correct reaction type given sources and sinks.
    If a molecule is provided, copies the information in sources and sinks
    to be attached to the new molecule.
    """

    #NOTE: sinks need to be used because there are dummy sources, but never dummy sinks.
    orig_mol =  sinks[0].molecule 
    if mol:
        new_sources = utils.shift_molecule_references(sources,mol)
        new_sinks = utils.shift_molecule_references(sinks,mol)
    else:
        new_sources = sources
        new_sinks = sinks
    if r_type == "proton_transfer":
        return ProtonTransfer(new_sources,new_sinks)
    elif r_type == "SN2":
        return SN2(new_sources,new_sinks)
    elif r_type == "AN":
        return AN(new_sources,new_sinks)
    elif r_type == "DN":
        return DN(new_sources,new_sinks)
    elif r_type == "E2":
        return E2(new_sources,new_sinks)
    elif r_type == "EB":
        return EB(new_sources,new_sinks)
    elif r_type == "E1":
        return E1(new_sources,new_sinks)
    elif r_type == "AE" and second_product == False:
        return AE(new_sources,new_sinks)
    elif r_type == "ADE3" and second_product == False:
        print("ADE3 1 was called ")
        return ADE3(new_sources,new_sinks)
    elif r_type == "ADN":
        return ADN(new_sources,new_sinks)
    elif r_type == "NuL":
        return NuL(new_sources,new_sinks)
    #If we need t oproduce the second_product for Elemination or Addition reactions if Mark rule is equal on both carbons
    elif r_type == "AE" and second_product == True:
        return AE(new_sources,new_sinks,second_product=True)
    elif r_type == "ADE3" and second_product == True:
        print("ADE3 2 was called ")
        return ADE3(new_sources,new_sinks,second_product=True)
       
    else:
        raise ValueError("Reaction type {} is not supported.".format(r_type))
