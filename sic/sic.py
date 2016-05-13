#!/usr/bin/python
#coding=utf-8
"""
This is the command-line interface for SiC³. It takes as input the same options as the
old SiC, and produces similar output on request.

This interface is provided purely for creating similar functionality to the old SiC
and inter-language cooperation; if other Python programs wish to use the functionality
of SiC³, they can import the libraries used here to produce results.

As of now, this program only calls the old SiC's command-line arguments.
"""
import argparse
import os #for making/destroying files for the sake of SiC, whose -I argument works a little oddly.
import re #for "arbitrary delimiter" support
import sys #for exit codes
from sic_io import sic_io#for parsing SiC-format input files
import logging #for debug logs - worry about this later
from brain import decision_engine

SIC_PATH = "/home/sic/sic/sic" #that's just sic.

def find_mechanism(reac,prod,solv=False):
    """
    Where the magic happens. Finds the mechanism by copying the current reaction state into a
    new set of Molecule objects, generating choices, and picking the best one.
    Most of the work is done outside of this module, but the core is left here so that
    other programs (such as SiGC) can access the full functionality without
    having to import a bunch of stuff.
    """
    #first convert them to our "state SMILES". SMILES is ASCII so make sure to change into that,
    #specifically because OBabel will choke on Unicode since it doesn't match std::string.
    reactants = sic_io.create_state_smiles(reac).encode("ascii","ignore")
    products = sic_io.create_state_smiles(prod).encode("ascii","ignore")
    solvent = sic_io.create_state_smiles(solv).encode("ascii","ignore") if solv else False
    try:
        mech = decision_engine.get_mechanism(reactants,products,solvent=solvent)
    except ValueError as e:
        print(e)
        return "Error encountered while writing up mechanism."
    return sic_io.write_up_mechanism(mech,solvent=solvent)

#not sure why you'd want to import this package, but it's good practice to wrap all argparse calls in this
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Command-line interface to SiC³")
    parser.add_argument("-i","--input-file",help="Path to the SiC-format input file. See SiC thesis for details. If \
            -r and -p are present, they will be prioritized over this argument.")
    parser.add_argument("-o","--output-file",help="Redirect all output produced by SiC³ to this file")
    parser.add_argument("-d","--debug",help="Include debug output")
    parser.add_argument("-r","--reactant",help="Reactant molecules, as a SMILES string. Takes in multiple arguments",action="append")
    parser.add_argument("-p","--product",help="Product molecules, as a SMILES string. Takes in multiple arguments",action="append")
    parser.add_argument("-s","--solvent",help="Solvent molecules, as a SMILES string. Takes in multiple arguments. \
            Not currently implemented, and will raise a NotImplementedError",action="append")
    parser.add_argument("-g","--graphics",help="Produces a graphical representation of reactant, product, solvent, and intermediate \
            molecules. Not currently implemented, and will raise a NotImplementedError",action="store_true")
    args = parser.parse_args()
    react_obj = False #will get filled in the if block below
    if args.solvent:
        raise NotImplementedError("Solvents aren't implemented yet.")
    if args.graphics:
        raise NotImplementedError("Graphical representations aren't implemented yet.")
    if args.reactant and args.product:
        react_obj = {"reactants": args.reactant, "products": args.product}
    else:
        if args.input_file:
            sic_input = open(args.input_file) #if there's an exception the user should see it, catching it does no good
            react_obj = sic_io.parse_sic_file(sic_input)
        else:
            print("Reactants and Products need to be provided, whether by input file or by arguments, in order for SiC³ to find a mechanism.")
            sys.exit(1)
    print find_mechanism(react_obj["reactants"],react_obj["products"],solv=(react_obj["solvent"] if react_obj.has_key("solvent") else False))



