#!/usr/bin/python
#coding=utf-8
"""
Tests whether halogen acids and their conjugate bases have their pKa correctly identified.
"""

import unittest
from .. import pka_chart

class HalogenPkaTest(unittest.TestCase):
	def setUp(self):
		self.PKA_CHART = pka_chart.PKA_CHART
		self.waters = ["[OH-]","O","[OH3+]"]
		self.fluorines = ["F","[F-]"]
		self.chlorines = ["Cl","[Cl-]"]
		self.bromines = ["Br","[Br-]"]
		self.iodines = ["I","[I-]"]
	
	def testWaterPka(self):


def main():
	unittest.main()

if __name__ == "__main__":
	main()
