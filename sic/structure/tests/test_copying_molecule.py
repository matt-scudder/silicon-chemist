"""
Tests whether we can make a copy of a molecule and do the following things:
    1. Make bond changes in the copy without affecting the original.
    2. Make bond changes using the indices in the original, and still modifying the correct atoms in the copy.
    3. Shift references on source/sink objects and produce the same indices on the other end.

It is crucially important that both of these work in order for our ReactionState to work properly.
"""

import unittest
import structure.struct_ops as struct_ops
import pybel
import segmentation.segmentation as segmentation
import pka.pka as pka
import utils

class MoleculeCopyAndModifyTest(unittest.TestCase):
    def setUp(self):
        prod_mol = pybel.readstring("smi","[F-].[OH3+]")
        prod_mol.addh()
        self.prod_mol = prod_mol
        orig_mol = pybel.readstring("smi","F.O")
        orig_mol.addh()
        pka.get_all_pka(orig_mol)
        segment = segmentation.segment_molecule(orig_mol)
        self.orig_mol = orig_mol
        self.sources = segment["sources"]
        self.sinks = segment["sinks"]
        #make mol copies in the tests involving modifying the copy, so we don't lose track of state


    def testShiftReference(self):
        """
        Tests whether utils.shift_molecule_references actually works.
        """
        copy_mol = struct_ops.copy_molecule(self.orig_mol)
        new_sources = utils.shift_molecule_references(self.sources,copy_mol)
        for source in new_sources:
            self.assertTrue(source.molecule == copy_mol)
            self.assertFalse(source.molecule == self.orig_mol)

    def testCopyNotModifiedAfterModifyingOriginal(self):
        """
        Tests the case of making a copy and modifying the original.
        Also, we regenerate the molecule here instead of using setUp because of
        state changes in orig_mol, since we don't know the order in which tests
        will be run.
        """
        orig_mol = pybel.readstring("smi","F.O")
        orig_mol.addh()
        pka.get_all_pka(orig_mol)
        #now get sources and sinks
        segment = segmentation.segment_molecule(orig_mol)
        #copy using struct_ops
        copy_mol = struct_ops.copy_molecule(orig_mol)
        source = segment["sources"][0]
        sink = segment["sinks"][0]
        struct_ops.make_bond(source.get_atom("Y"),sink.get_atom("H"),orig_mol)
        struct_ops.break_bond(sink.get_atom("H"),sink.get_atom("L"),orig_mol)
        self.assertTrue(orig_mol.write("can") == self.prod_mol.write("can"))
        self.assertFalse(orig_mol.write("can") == copy_mol.write("can"))

    def testOriginalNotModifiedAfterModifyingCopy(self):
        """
        Tests the far more usual case of making the copy and then modifying it.
        This requires that indices are preserved when making the copy.
        """
        #first make the copy
        copy_mol = struct_ops.copy_molecule(self.orig_mol)
        #shift references using utility function
        new_sources = utils.shift_molecule_references(self.sources,copy_mol)
        new_sinks = utils.shift_molecule_references(self.sinks,copy_mol)
        Y = new_sources[0].get_atom("Y")
        H = new_sinks[0].get_atom("H")
        L = new_sinks[0].get_atom("L")
        struct_ops.make_bond(Y,H,copy_mol)
        struct_ops.break_bond(H,L,copy_mol)
        self.assertTrue(copy_mol.write("can") == self.prod_mol.write("can"))
        self.assertFalse(self.orig_mol.write("can") == copy_mol.write("can"))

if __name__ == "__main__":
    unittest.main()
