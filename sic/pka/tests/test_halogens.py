#!/usr/bin/python
#coding=utf-8
"""
Tests whether halogen acids and their conjugate bases have their pKa correctly identified.
"""

import unittest
from .. import pka_chart,pka
from pybel import readstring

class HalogenPkaTest(unittest.TestCase):
	def setUp(self):
		self.PKA_CHART = pka_chart.PKA_CHART
		self.fluorines = ["F","[F-]"]
		self.chlorines = ["Cl","[Cl-]"]
		self.bromines = ["Br","[Br-]"]
		self.iodines = ["I","[I-]"]
	
	def testFluorinePka(self):
		fluorine_mols = []
		for s_string in self.fluorines:
			mol = readstring("smi",s_string)
			mol.addh()
			pka.get_all_pka(mol)
			fluorine_mols.append(mol)
		for atom in fluorine_mols[0]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,fluorine_mols[0]) == 3.2)
			else:
				self.assertTrue(pka.get_pka(atom.idx,fluorine_mols[0]) == None)
		for atom in fluorine_mols[1]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,fluorine_mols[1]) == None)
			else:
				self.assertTrue(pka.get_pka(atom.idx,fluorine_mols[1]) == 3.2)

	def testChlorinePka(self):
		chlorine_mols = []
		for s_string in self.chlorines:
			mol = readstring("smi",s_string)
			mol.addh()
			pka.get_all_pka(mol)
			chlorine_mols.append(mol)
		for atom in chlorine_mols[0]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,chlorine_mols[0]) == -7)
			else:
				self.assertTrue(pka.get_pka(atom.idx,chlorine_mols[0]) == None)
		for atom in chlorine_mols[1]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,chlorine_mols[1]) == None)
			else:
				self.assertTrue(pka.get_pka(atom.idx,chlorine_mols[1]) == -7)

	def testBrominePka(self):
		bromine_mols = []
		for s_string in self.bromines:
			mol = readstring("smi",s_string)
			mol.addh()
			pka.get_all_pka(mol)
			bromine_mols.append(mol)
		for atom in bromine_mols[0]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,bromine_mols[0]) == -9)
			else:
				self.assertTrue(pka.get_pka(atom.idx,bromine_mols[0]) == None)
		for atom in bromine_mols[1]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,bromine_mols[1]) == None)
			else:
				self.assertTrue(pka.get_pka(atom.idx,bromine_mols[1]) == -9)

	def testIodinePka(self):
		iodine_mols = []
		for s_string in self.iodines:
			mol = readstring("smi",s_string)
			mol.addh()
			pka.get_all_pka(mol)
			iodine_mols.append(mol)
		for atom in iodine_mols[0]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,iodine_mols[0]) == -10)
			else:
				self.assertTrue(pka.get_pka(atom.idx,iodine_mols[0]) == None)
		for atom in iodine_mols[1]:
			if atom.atomicnum == 1:
				self.assertTrue(pka.get_pka(atom.idx,iodine_mols[1]) == None)
			else:
				self.assertTrue(pka.get_pka(atom.idx,iodine_mols[1]) == -10)




def main():
	unittest.main()

if __name__ == "__main__":
	main()
