/*
 * Copyright (C) 2003-2015 Syed Asad Rahman <asad @ ebi.ac.uk>.
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
package uk.ac.ebi.reactionblast.mapping.algorithm;

import static java.lang.System.out;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.ISAROMATIC;
import static org.openscience.cdk.CDKConstants.ISINRING;
import static org.openscience.cdk.CDKConstants.RING_CONNECTIONS;
import static org.openscience.cdk.DefaultChemObjectBuilder.getInstance;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import static org.openscience.cdk.graph.Cycles.all;
import static org.openscience.cdk.graph.Cycles.or;
import static org.openscience.cdk.graph.Cycles.relevant;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import static org.openscience.cdk.interfaces.IBond.Order.SINGLE;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.smsd.helper.MoleculeInitializer;

/**
 *
 * @author Asad
 */
class IsomeraseHandler {

    private static final Logger LOG = getLogger(IsomeraseHandler.class.getName());

    /**
     *
     * @param atomContainer Atom container where rings are to be marked
     * @throws CDKException if there is a problem in ring perception or
     * aromaticity detection, which is usually related to a timeout in the ring
     * finding code.
     */
    protected static void initializeMolecule(IAtomContainer atomContainer) throws CDKException {
        MoleculeInitializer.initializeMolecule(atomContainer);
    }

    private final Map<IRingSet, IAtomContainer> ringContainerCountR;
    private final Map<IRingSet, IAtomContainer> ringContainerCountP;
    private final IRingSet sssrEduct;
    private final IRingSet sssrProduct;
    private final boolean DEBUG = false;

    protected final IReaction reaction;

    IsomeraseHandler(IReaction reaction) {
        this.reaction = reaction;
        ringContainerCountR = getRingContainerCount(reaction.getReactants());
        ringContainerCountP = getRingContainerCount(reaction.getProducts());

        if (ringContainerCountR.size() == 1 && ringContainerCountP.size() == 1) {
            IAtomContainer educt = ringContainerCountR.values().iterator().next();
            IAtomContainer product = ringContainerCountP.values().iterator().next();
            try {
                initializeMolecule(educt);
            } catch (CDKException ex) {
                getLogger(GameTheoryRings.class.getName()).log(SEVERE, null, ex);
            }
            try {
                initializeMolecule(product);
            } catch (CDKException ex) {
                getLogger(GameTheoryRings.class.getName()).log(SEVERE, null, ex);
            }
            // sets SSSR information
            SSSRFinder finder = new SSSRFinder(educt);
            sssrEduct = finder.findEssentialRings();
            finder = new SSSRFinder(product);
            sssrProduct = finder.findEssentialRings();
            boolean chipBondInTheRing = chipBondInTheRing(educt, product);

        } else {
            sssrEduct = null;
            sssrProduct = null;
            IAtomContainerSet reactants = reaction.getReactants();
            IAtomContainerSet products = reaction.getProducts();

            if (reactants.getAtomContainerCount() == 1 && products.getAtomContainerCount() == 1) {
                try {
                    IAtomContainer educt = reactants.atomContainers().iterator().next();
                    IAtomContainer product = products.atomContainers().iterator().next();
                    boolean chipPhophateInSingleReactantProductNotInRing = chipPhophateInSingleReactantProductNotInRing(educt, product);
                } catch (CDKException ex) {
                    getLogger(IsomeraseHandler.class.getName()).log(SEVERE, null, ex);
                }
            }
        }
    }

    /*
     * Returns Number of Container with Ring Atoms
     */
    private Map<IRingSet, IAtomContainer> getRingContainerCount(IAtomContainerSet acs) {
        CycleFinder cycles = or(all(), all());
        Map<IRingSet, IAtomContainer> ringSet = new HashMap<>();
        for (IAtomContainer ac : acs.atomContainers()) {
            IRingSet basicRings;
            try {
                Cycles rings = cycles.find(ac);
                basicRings = rings.toRingSet();
                if (!basicRings.isEmpty()) {
                    basicRings.setID(ac.getID());
                    ringSet.put(basicRings, ac);
                }
            } catch (Intractable ex) {
                getLogger(IsomeraseHandler.class.getName()).log(SEVERE, null, ex);
            }
        }
        return ringSet;
    }

    private boolean chipBondInTheRing(IAtomContainer educt, IAtomContainer product) {

        if (sssrEduct.getAtomContainerCount() == 1
                && sssrProduct.getAtomContainerCount() == 1) {
            IAtomContainer ringE = sssrEduct.atomContainers().iterator().next();
            IAtomContainer ringP = sssrProduct.atomContainers().iterator().next();

            if (ringE.getAtomCount() == 6 && ringP.getAtomCount() == 5) {
                return findAndChipBond(ringE, educt);
            } else if (ringE.getAtomCount() == 5 && ringP.getAtomCount() == 6) {
                return findAndChipBond(ringP, product);
            }
        } else if (sssrEduct.getAtomContainerCount() == 2
                && sssrEduct.getAtomContainerCount() == sssrProduct.getAtomContainerCount()) {
            return findAndChipBondBetweenRings(educt) | findAndChipBondBetweenRings(product);
        }

        return false;
    }

    /*
     * Example reaction is R01518
     */
    private boolean chipPhophateInSingleReactantProductNotInRing(IAtomContainer educt, IAtomContainer product) throws CDKException {

        if (ringContainerCountR.isEmpty() && ringContainerCountP.isEmpty()) {
            String phosphateSMILES = "OP(O)(O)=O";
            SMARTSQueryTool smartsPhosphate;

            if (DEBUG) {
                out.println("String phosphateSMILES = \"OP(O)(O)=O\";");
            }

            smartsPhosphate = new SMARTSQueryTool(phosphateSMILES, getInstance());
            boolean matchesE = smartsPhosphate.matches(educt);
            boolean matchesP = smartsPhosphate.matches(product);

            if (DEBUG) {
                out.println("String smartsPhosphate Q " + matchesE);
                out.println("String smartsPhosphate T " + matchesP);
            }
            if (matchesE && matchesP) {
                boolean findAndChipBondPhophate1 = findAndChipBondPhophate(educt);
                boolean findAndChipBondPhophate2 = findAndChipBondPhophate(product);
                if (DEBUG) {
                    out.println("findAndChipBondPhophate1 Q " + findAndChipBondPhophate1);
                    out.println("findAndChipBondPhophate2 T " + findAndChipBondPhophate2);
                    out.println("SM " + new SmilesGenerator().create(educt));
                    out.println("SM " + new SmilesGenerator().create(product));
                }
                return findAndChipBondPhophate1 && findAndChipBondPhophate2;
            }
        }

        return false;
    }

    private boolean findAndChipBond(IAtomContainer container, IAtomContainer referenceContainer) {
        boolean flag = false;

        if (container != null) {
            for (IBond bond : container.bonds()) {
                if ((bond.getAtom(0).getSymbol().equalsIgnoreCase("O")
                        && bond.getAtom(1).getSymbol().equalsIgnoreCase("C"))
                        || (bond.getAtom(0).getSymbol().equalsIgnoreCase("C")
                        && bond.getAtom(1).getSymbol().equalsIgnoreCase("O"))) {
                    if (!bond.getAtom(0).getFlag(ISAROMATIC)
                            && !bond.getAtom(1).getFlag(ISAROMATIC)) {
                        if (referenceContainer.contains(bond)) {
                            IAtom atom = bond.getAtom(0).getSymbol().equalsIgnoreCase("C") ? bond.getAtom(0) : bond.getAtom(1);
                            List<IBond> neighbourhoodBonds = referenceContainer.getConnectedBondsList(atom);
                            flag = false;
                            for (IBond neighbourhoodBond : neighbourhoodBonds) {
                                if (neighbourhoodBond.contains(atom) && !neighbourhoodBond.getFlag(ISINRING)) {
                                    if ((neighbourhoodBond.getAtom(0).getSymbol().equalsIgnoreCase("O")
                                            && neighbourhoodBond.getAtom(1).getSymbol().equalsIgnoreCase("C"))
                                            || (neighbourhoodBond.getAtom(0).getSymbol().equalsIgnoreCase("C")
                                            && neighbourhoodBond.getAtom(1).getSymbol().equalsIgnoreCase("O"))) {
                                        flag = true;
                                    }
                                }
                            }

                            if (flag) {
                                referenceContainer.removeBond(bond);
                                break;
                            }
                        }
                    }
                }
            }
        }
        return flag;
    }

    private boolean findAndChipBondPhophate(IAtomContainer container) {
        boolean flag = false;

        for (IBond bond : container.bonds()) {
            IAtom atomE = bond.getAtom(0);
            IAtom atomP = bond.getAtom(1);
            if ((atomE.getSymbol().equals("O") && atomP.getSymbol().equals("P"))
                    || (atomE.getSymbol().equals("P") && atomP.getSymbol().equals("O"))) {
                if (bond.getOrder().equals(SINGLE)) {
                    IAtom oxygen = atomE.getSymbol().equals("O") ? atomE : atomP;
                    List<IBond> neighbourBonds = container.getConnectedBondsList(oxygen);
                    if (neighbourBonds.size() == 2) {
                        for (IBond b : neighbourBonds) {
                            if (b.getAtom(0).getSymbol().equals("O")
                                    || b.getAtom(0).getSymbol().equals("P")) {
                                if (b.getAtom(1).getSymbol().equals("O")
                                        || b.getAtom(1).getSymbol().equals("P")) {
                                    container.removeBond(b);

                                    if (DEBUG) {
                                        out.println("bondToBeChipped " + b.getAtom(0).getSymbol());
                                        out.println("bondToBeChipped " + b.getAtom(1).getSymbol());
                                        out.println("removeBond o-p ");
                                    }
                                }
                            }
                        }
                        return true;
                    }
                }
            }
        }

        return flag;
    }

    /*
    This Method will find and chip the bonds between the rings
    Example R01557
     */
    private boolean findAndChipBondBetweenRings(IAtomContainer container) {
        if (DEBUG) {
            out.println("Find and Chip Bond Between Rings");
        }
        List<IBond> bond_to_be_removed = new ArrayList<>();
        for (IAtom atom : container.atoms()) {
            if (atom.getSymbol().equals("O")) {
                int number_of_rings = ((Integer) atom.getProperty(RING_CONNECTIONS));
                if (DEBUG) {
                    out.println("number_of_rings " + number_of_rings);
                }

                List<IBond> bonds = container.getConnectedBondsList(atom);
                if (DEBUG) {
                    out.println("number_of_bonds " + bonds.size());
                    for (IBond bond : bonds) {
                        out.println("BONDS "
                                + " B0 " + bond.getAtom(0).getSymbol()
                                + " B1 " + bond.getAtom(1).getSymbol());
                    }
                }

                if (bonds.size() == 2 && number_of_rings == 2) {
                    IBond bondToBeChipped = bonds.iterator().next();
                    bond_to_be_removed.add(bondToBeChipped);
                }
            }
        }
        for (IBond bond : bond_to_be_removed) {
            container.removeBond(bond);
            if (DEBUG) {
                try {
                    out.println("CHIPPING BONDS "
                            + " B0 " + bond.getAtom(0).getSymbol()
                            + " B1 " + bond.getAtom(1).getSymbol());
                    out.println("CHIPPED SM " + new SmilesGenerator().create(container));
                } catch (CDKException ex) {
                    getLogger(IsomeraseHandler.class.getName()).log(SEVERE, null, ex);
                }
            }
        }
        return !bond_to_be_removed.isEmpty();
    }

}
