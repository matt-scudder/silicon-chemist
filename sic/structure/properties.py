"""
Methods for getting properties of a structure rather than performing operations on it.
Examines the degree of a carbon, for example.
"""
from subprocess import run
from pathlib import Path
import utils
import re

# FileName of the Reaction Decoder JAR, located next to runserver.py, up one directory from here.
REACTION_DECODER_JAR = "rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar"

def get_mapping(reactants,products):
    """
    Gets the mapping between reactants and products.
    TODO: Replace this with an algorithm that doesn't rely on ReactionDecoder.
    Returns also the "ReactionDecoder canonical form" of reactants and products, since OpenBabel's canonical SMILES code
    fails for charged species.
    """
    mapping = {}
    #write out reactants/products string
    input_smiles = "%s>>%s" % (utils.write_mol(reactants),utils.write_mol(products))
    jar_path = f"{Path(__file__).parent.parent}/{REACTION_DECODER_JAR}"
    if not Path(jar_path).exists():
        raise RuntimeError("Can not find ReactionDecoder JAR", jar_path)
    rdt_process = run(["java","-jar",jar_path,"-Q","SMI","-q",input_smiles,"-g","-j","AAM","-f","TEXT"], capture_output=True, text=True)
    mapping_string = None
    if rdt_process.stderr:
        raise RuntimeError("Error from Reaction Decoder", rdt_process.stderr)
    if rdt_process.stdout:
        path_prefix = "Output is presented in text format: "
        AAM_text_file_path = None
        for stdout_line in rdt_process.stdout.split('\n'):
            if stdout_line.startswith(path_prefix):
                AAM_text_file_path = stdout_line[len(path_prefix):]
                break
        # get mapping string out of ECBLAST_smiles_AAM.txt file
        if not AAM_text_file_path:
            message = f"No output text file path found in RDT stdout. Looking for \"{path_prefix}\""
            raise RuntimeError(message, rdt_process.stdout)
        with open(AAM_text_file_path) as f:
            file_content = f.readlines()
            line_before_aam_mapping = "SELECTED AAM MAPPING\n"
            mapping_index = file_content.index(line_before_aam_mapping) + 1
            mapping_string = file_content[mapping_index]
            print(line_before_aam_mapping + mapping_string)
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
        for i in range(len(react_groups)):
            current_group = react_groups[i]
            number = int(current_group.split(":")[-1].replace("]","")) 
            internal_mapping[number] = i+1 #this should usually be the same...
        #now map against products
        for j in range(len(prod_groups)):  
            current_group = prod_groups[j]
            number = int(current_group.split(":")[-1].replace("]",""))
            if number <= (i+1) :  # check if the number of atoms in the reactants is the same in the products
                mapping[internal_mapping[number]] = j+1 #j+1 is the atom index in the product. internal_mapping[number] is the atom index in the reactant
            else:
                print("Oops! Wrong input. Please try again -_-")
                break
    else:
        raise ValueError("Could not find map between reactants and products.")
    print("Mapping is made ", (mapping))
    return (mapping,react_map,prod_map) #so that we can use the actual strings as our input and avoid issues with OBabel's buggy SMILES code


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
        if not (mol.OBMol.GetAtom(bond).IsHydrogen()):
            #H bonds don't count, neither do L if any
            #Update the above conditional for species other than L that don't count
            #though I'm fairly sure there aren't any.
            carbon_count += 1
    return carbon_count

def get_H_bonds(atom,mol):
    """
    Takes in an atom index and a Molecule object, and gets the number of hydrogen atoms bonded to it.
    You can use len() of this to get the number of hydrogens attached to an atom, a sort of
    reverse of get_carbon_degree.
    Intended to be used as a utility function for eliminations, and get_adjacent_ch.
    Assumes carbon, but if atom_label is set, can check for all H bonded to any provided generic
    class label.
    """
    H_bonds = []
    for bonded_atom in mol.connectivity_table.get_atoms_bonded(atom):
        if mol.OBMol.GetAtom(bonded_atom).IsHydrogen():
            H_bonds.append(bonded_atom)
    return H_bonds

def get_adjacent_ch(s_obj,carbon_label=False):
    """
    Takes in a source/sink object and returns the atom index of the carbon of an adjacent CH, if any,
    or returns False.
    Used to determine whether eliminations are possible.
    If carbon_label is set, uses that to search for the carbon instead of the string "C".
    """
    #NOTE: Possibly upgrade this to handle NH, OH, and so on?
    #NOTE: This can be made faster by not calling get_H_bonds, but it's fancier and better designed this way.
    valid_atomicnums = set([6])
    carb_string = carbon_label if carbon_label else "C"
    mol = s_obj.molecule
    obmol = mol.OBMol #typing is hard
    bonded_atom_dic = {} # a dictionary that stores the adjacent carbon index as a key, and the number of bonded H atoms to this carbon as a value
    for bonded_atom in mol.connectivity_table.get_atoms_bonded(s_obj.get_atom(carb_string)):
        atom = obmol.GetAtom(bonded_atom)
        if atom.GetAtomicNum() == 6: # checks if the adjacent atom is carbon
            if len(get_H_bonds(bonded_atom,mol)) > 0:
                return bonded_atom #return first one found
    return False 

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
        result_table[frozenset(atomlist)] = table[bond] #TODO: bug? what if the two atoms don't have a bond in the product?
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
    #print "mapping =", mapping
    #print "After remap: {}".format(mol1_bonds)
    #print "Product =", mol2_bonds
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

