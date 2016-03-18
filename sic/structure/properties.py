#!/usr/bin/python
#coding=utf-8
"""
Methods for getting properties of a structure rather than performing operations on it.
Examines the degree of a carbon, for example.
"""
import subprocess
import utils
import re

REACTION_DECODER_PATH = "/home/sic/newsic/sic/ReactionDecoder/bin" #TODO: Figure out how to not need this

def get_mapping(reactants,products):
    """
    Gets the mapping between reactants and products.
    TODO: Replace this with an algorithm that doesn't rely on ReactionDecoder.
    """
    mapping = {}
    #write out reactants/products string
    input_smiles = "%s>>%s" % (utils.write_mol(reactants),utils.write_mol(products))
    AAM_output = subprocess.check_output(["java","-jar","%s/usableRDT.jar"%REACTION_DECODER_PATH,"-Q","SMI","-q",input_smiles,"-g","-j","AAM","-f","TEXT"])
    found_mapping = False
    mapping_string = None
    for line in AAM_output.split("\n"):
        if "SELECTED AAM MAPPING" in line:
            found_mapping = True
            continue
        if found_mapping:
            mapping_string = line
            break
    if mapping_string:
        #remove hydrogen-only maps - these are unlikely but they do come up and mess up the rest of our procedure
        #also I'm not digging into their code to figure out why they do this only *sometimes*, I'm mad
        #enough at Java as it is.
        hydrogen_re = re.compile("\[H:[0-9]*\]")
        mapping_string = hydrogen_re.sub("",mapping_string)
        #split into reactant and product map
        react_map,prod_map = mapping_string.split(">>")
        #look for the [] groups
        react_groups = re.findall("\[[^\]]+\]",react_map)
        prod_groups = re.findall("\[[^\]]+\]",prod_map)
        internal_mapping = {} #to make up for the H groups we excised. Most of the time this will just be a map from a number to itself, but gotta future-proof.
        for i in xrange(len(react_groups)):
            current_group = react_groups[i]
            number = int(current_group.split(":")[-1].replace("]",""))
            internal_mapping[number] = i+1 #this should usually be the same...
        #now map against products
        for j in xrange(len(prod_groups)):
            current_group = prod_groups[j]
            number = int(current_group.split(":")[-1].replace("]",""))
            mapping[internal_mapping[number]] = j+1 #j+1 is the atom index in the product. internal_mapping[number] is the atom index in the reactant
    return mapping 


def get_bonds(atom,mol):
    """
    Takes in an atom index and its matching molecule, and gets the indices of the atoms that it's bound to, based on the connectivity table.
    Utility function for less typing and clearer meanings.
    """
    return mol.connectivity_table[atom]


def get_carbon_degree(s_obj,carbon_label=False):
    """
    Takes in an source/sink object, and gets the number of non-H atoms attached to its carbon.
    Used to determine whether a carbon is primary, secondary, tertiary or methyl.
    If carbon_label is set, uses that to search for the carbon instead of the string "C".
    """
    carb_string = carbon_label if carbon_label else "C"
    carbon_count = 0
    mol = s_obj.molecule
    has_L = "L" in s_obj.atoms
    for bond in mol.connectivity_table.get_atoms_bonded(s_obj.get_atom(carb_string)):
        if not (mol.OBMol.GetAtom(bond).IsHydrogen()) and (has_L and bond != s_obj.get_atom("L")):
            #H bonds don't count, neither do L if any
            #Update the above conditional for species other than L that don't count
            #though I'm fairly sure there aren't any.
            carbon_count += 1
    return carbon_count

def remap_bonds(table,mapping):
    """
    Takes as input a closer_to_product_table, and returns a set of frozenset bonds
    that have been mapped.
    """
    result_table = {}
    for bond in table:
        atomlist = []
        for atom in bond:
            if atom == "H":
                atomlist.append(atom)
            else:
                atomlist.append(mapping[atom])
        result_table[frozenset(atomlist)] = table[bond]
    return result_table

def get_bond_distance(mol1,mol2,mapping):
    """
    Gets the difference in bonds between two Molecule objects, using a 1:1 atom-to-atom mapping between them.
    This mapping is assumed to be a dictionary with the following key:value format:

    {atom_idx_in_mol1 : atom_idx_in_mol2}

    Bond distance is calculated by checking how many match. If the distance is lower (closer to 0), then mol1 is close to mol2.
    If this function returns 0, mol1 == mol2, chemically.
    
    This function assumes both Molecule objects have had a connectivity_table generated for them (see generate_connectivity_table above).
    """
    #TODO: Replace all prints with logging.debug stuff
    difference_counter = 0
    #get all bonds in mol1 and in mol2 using our closer_to_product table
    #print "Before remap: {}".format(mol1.connectivity_table.closer_to_product_table)
    mol1_bonds = remap_bonds(mol1.connectivity_table.closer_to_product_table,mapping)
    mol2_bonds = mol2.connectivity_table.closer_to_product_table
    #print "After remap: {}".format(mol1_bonds)
    #print "Product bonds: {}".format(mol2_bonds)
    for bond in mol1_bonds:
        if bond not in mol2_bonds:
            #if the bond is not present at all, add the bond order of these non-present bonds
            #print "Bond not present: %s" %bond
            difference_counter += mol1_bonds[bond]
        else:
            diff = abs(mol2_bonds[bond] - mol1_bonds[bond])
            #if the bond is in fact present, add the absolute bond difference, as it is possible that mol2
            #print "Bond is present: %s, difference is %s" % (bond,diff)
            #has less bonds than mol1 at a particular center (e.g. additions).
            difference_counter += diff
    #now check the other direction
    for bond in mol2_bonds:
        if bond not in mol1_bonds:
            #if the bond is present in mol2 but not in mol1, it's also a difference.
            #print "Bond not present: %s" % bond
            difference_counter += mol2_bonds[bond] 
    return difference_counter
