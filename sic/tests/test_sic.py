from .. import sic
import unittest

PARSE_FILE_PATHS = ["reaction_files/in.1", "reaction_files/in.9"]
SOLVENT_FILE = "reaction_files/in.8"

class SiCParserTest(unittest.TestCase):
	def setUp(self):
		self.commafile = open(PARSE_FILE_PATHS[0],"r")
		self.periodfile = open(PARSE_FILE_PATHS[1],"r")
		self.solvfile = open(SOLVENT_FILE,"r")

	def testCommaParse(self):
		"""
		Tests whether we parse comma-delimited files correctly.
		"""
		first_reaction = sic.parse_sic_file(self.commafile)
		self.assertEqual(first_reaction["reactants"],["[OH3+]","[OH3+]","[O-]C(=O)CC([O-])CCC(=O)[O-]"])
		self.assertEqual(first_reaction["products"],["O","O",'OC(=O)CC(OH)CCC(=O)[O-]'])
		self.assertFalse(first_reaction["solvent"]) #there shouldn't be solvent here

	def testPeriodParse(self):
		"""
		Tests whether we parse period-delimited files correctly.
		"""
		first_reaction = sic.parse_sic_file(self.periodfile)
		self.assertEqual(first_reaction["reactants"],["C=CC(C)(O)CCC=C(C)C","I"])
		self.assertEqual(first_reaction["products"],["CC1=CCC(C(C)=C)CC1","[OH3+]","[I-]"])
		self.assertFalse(first_reaction["solvent"]) #there shouldn't be solvent here

	def testSolventParse(self):
		"""
		Tests whether we read in solvent. Also tests mixed use of delimiters.
		"""
		first_reaction = sic.parse_sic_file(self.solvfile)
		self.assertEqual(first_reaction["reactants"],["C1CCC(O)C(C1)(C2)CCC2"])
		self.assertEqual(first_reaction["products"],["C1CCCC(=C12)CCCC2","O","[O-]S(=O)(=O)[O-]"])
		self.assertEqual(first_reaction["solvent"],["OS(=O)(=O)O","O"]) 

def main():
	unittest.main()

if __name__ == "__main__":
	main()
