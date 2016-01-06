"""
This is the command-line interface for SiC³. It takes as input the same options as the
old SiC, and produces similar output on request.

This interface is provided purely for creating similar functionality to the old SiC
and inter-language cooperation; if other Python programs wish to use the functionality
of SiC³, they can import the libraries used here to produce results.

As of now, this program only calls the old SiC's command-line arguments.
"""

SIC_PATH = "/home/sic/sic/sic" #that's just sic.

import argparse
import os

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


def parse_sic_file(sic_input):
	"""
	Parses SiC-format file and returns an object with reactants, solvent and products.

	SiC format is defined by the following pattern:

	{reactants}>>{solvent}>>{products}

	With delimiters between each reactant or product being arbitrary, SO LONG AS
	the same delimiter is used between each reactant, product or solvent.
	Characters that CANNOT be used as delimiters are SMILES characters such as =, #,
	[,], etc. (see SMILES format guide at http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html ).

	This function Separates out reactants by taking in all the characters before the first >,
	solvent by taking in any characters found after the first two > in the string,
	and products by taking in all the characters after the last >.
	"""
	return {}
