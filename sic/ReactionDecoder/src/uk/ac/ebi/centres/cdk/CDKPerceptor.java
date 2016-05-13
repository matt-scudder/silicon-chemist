/*
 * Copyright (c) 2012. John May
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */
package uk.ac.ebi.centres.cdk;

import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.exception.CDKException;
import static org.openscience.cdk.geometry.GeometryTools.has2DCoordinates;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import uk.ac.ebi.centres.DefaultPerceptor;
import uk.ac.ebi.centres.priority.AtomicNumberRule;
import uk.ac.ebi.centres.priority.CombinedRule;
import uk.ac.ebi.centres.priority.DuplicateAtomRule;
import uk.ac.ebi.centres.priority.MassNumberRule;
import uk.ac.ebi.centres.priority.access.AtomicNumberAccessor;
import uk.ac.ebi.centres.priority.access.MassNumberAccessor;
import uk.ac.ebi.centres.priority.access.PsuedoAtomicNumberModifier;
import uk.ac.ebi.centres.priority.access.descriptor.AuxiliaryDescriptor;
import uk.ac.ebi.centres.priority.access.descriptor.PrimaryDescriptor;
import uk.ac.ebi.centres.priority.descriptor.PairRule;
import uk.ac.ebi.centres.priority.descriptor.RSRule;
import uk.ac.ebi.centres.priority.descriptor.ZERule;

/**
 * @author John May
 */
public class CDKPerceptor extends DefaultPerceptor<IAtom> {

    private final static ILoggingTool logger
            = createLoggingTool(CDKPerceptor.class);

    private static final Logger LOG = getLogger(CDKPerceptor.class.getName());

    /**
     *
     */
    public CDKPerceptor() {
        super(new CombinedRule<>(
                new AtomicNumberRule<>(
                        new PsuedoAtomicNumberModifier<>(
                                new AtomicNumberAccessor<IAtom>() {

                                    @Override
                                    public int getAtomicNumber(IAtom atom) {
                                        /*
                                        * if its null, assuming its an "R" then put the masss less than element Carbon (6) TO DO this
                                        * fix properly
                                        */
                                        return atom.getAtomicNumber() == null ? 0 : atom.getAtomicNumber();
                                    }
                                })),
                new DuplicateAtomRule<IAtom>(),
                new MassNumberRule<>(new MassNumberAccessor<IAtom>() {

                    @Override
                    public int getMassNumber(IAtom atom) {
                        /*
                        * if its null, assuming its an "R" then put the masss less than element Carbon
                        */
                        return atom.getMassNumber() == null ? 11 : atom.getMassNumber();
                    }
                }),
                new ZERule<IAtom>(),
                new PairRule<>(new PrimaryDescriptor<IAtom>()),
                new RSRule<>(new PrimaryDescriptor<IAtom>())),
                new CombinedRule<>(
                        new AtomicNumberRule<>(
                                new PsuedoAtomicNumberModifier<>(
                                        new AtomicNumberAccessor<IAtom>() {

                                            @Override
                                            public int getAtomicNumber(IAtom atom) {
                                                /*
                                                * if its null, assuming its an "R" then put the masss less than element Carbon TO DO this fix
                                                * properly
                                                */
                                                return atom.getAtomicNumber() == null ? 0 : atom.getAtomicNumber();
                                            }
                                        })),
                        new MassNumberRule<>(new MassNumberAccessor<IAtom>() {

                            @Override
                            public int getMassNumber(IAtom atom) {
                                return atom.getMassNumber();
                            }
                        }),
                        new ZERule<IAtom>(),
                        new PairRule<>(new AuxiliaryDescriptor<IAtom>()),
                        new RSRule<>(new AuxiliaryDescriptor<IAtom>())),
                new CDK2DSignCalculator());
    }

    /**
     *
     * @param container
     */
    public void perceive(IAtomContainer container) {
        try {
            /*
            Check for 2D co-ordinates for EC-BLAST, must else it will fail!
            */
            if (!has2DCoordinates(container)) {
                try {
                    /*
                    Clone it else it will loose mol ID
                    */
                    IAtomContainer clone = container.clone();
                    StructureDiagramGenerator sdg = new StructureDiagramGenerator(clone);
                    sdg.generateCoordinates();
                    container = sdg.getMolecule();
                } catch (CDKException e) {
                }
            }
            perceive(new CDKCentreProvider(container), new CDKManager(container));
        } catch (Exception e) {
            logger.warn("WARNING: 2D CDK based stereo perception failed! ");
        }
    }
}
