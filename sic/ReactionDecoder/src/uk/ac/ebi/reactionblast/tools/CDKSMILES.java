/*
 * Copyright (C) 2007-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.tools;

import java.io.IOException;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import static org.openscience.cdk.smiles.CanonSmiAdapter.create;
import org.openscience.cdk.smiles.SmilesGenerator;
import static org.openscience.cdk.smiles.SmilesGenerator.generic;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeDayLight;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.cloneWithIDs;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 *
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class CDKSMILES {
    private static final Logger LOG = getLogger(CDKSMILES.class.getName());

    private final ILoggingTool logger
            = createLoggingTool(CDKSMILES.class);
    private final IAtomContainer molecule;

    /**
     *
     * @param mol
     * @param removeH
     * @param remove_AAM
     * @throws CloneNotSupportedException
     */
    public CDKSMILES(IAtomContainer mol, boolean removeH, boolean remove_AAM) throws CloneNotSupportedException {
        this.molecule = removeHydrogensExceptSingleAndPreserveAtomID(cloneWithIDs(mol));
        if (remove_AAM) {
            for (IAtom a : molecule.atoms()) {
                a.removeProperty(ATOM_ATOM_MAPPING);
            }
        }
    }

    /**
     *
     * @return
     */
    public String getCanonicalSMILES() {
        String smiles = "NA";
        if (molecule.getAtomCount() == 0) {
            return smiles;
        }
        try {
            aromatizeDayLight(molecule);
            return create(molecule);
        } catch (CDKException ex) {
            logger.error("ERROR : in generating CDK SMILES" + molecule.getID());
        } catch (IOException ex) {
            getLogger(CDKSMILES.class.getName()).log(SEVERE, null, ex);
        }
//        /*
//         Chemaxon based unique SMILES generator for debugging
//         */
//        try {
//            chemaxon.struc.Molecule chemAxonMolecule = uk.ac.ebi.tools.chemaxonhandler.CDKChemaxonIOConveter.getChemAxonMolecule(molecule);
//            /*
//             * Set Absolute stereo function
//             */
//            chemAxonMolecule.setAbsStereo(true);
//
//            /*
//             * check valences must before hydrogenize to calculate implicit H
//             *
//             */
//            chemAxonMolecule.valenceCheck();
//
//            /*
//             * check Hybridization must before hydrogenize to calculate implicit H
//             *
//             */
//            chemAxonMolecule.calcHybridization();
//            /*
//             * check aromaticity must before hydrogenize to calculate implicit H
//             *
//             */
//            chemAxonMolecule.aromatize();
//            return chemaxon.formats.MolExporter.exportToFormat(chemAxonMolecule, "smiles:au");
//        } catch (Exception ex) {
//            Logger.getLogger(CDKSMILES.class.getName()).log(Level.SEVERE, null, ex);
//        }

        return smiles;
    }

    /**
     *
     * @return
     */
    public String getGenericSMILES() {
        String smiles = "NA";
        if (molecule.getAtomCount() == 0) {
            return smiles;
        }
        try {
            SmilesGenerator g;
            g = generic().aromatic();
            return g.create(molecule);
        } catch (CDKException ex) {
            logger.error("ERROR : in generating CDK SMILES" + molecule.getID());
        }

        return smiles;
    }

    private boolean isPseudoAtoms() {
        for (IAtom atoms : molecule.atoms()) {
            if (atoms instanceof IPseudoAtom || atoms instanceof PseudoAtom) {
                return true;
            }
        }
        return false;
    }
}
