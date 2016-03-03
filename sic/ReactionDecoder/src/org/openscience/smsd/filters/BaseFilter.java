/* Copyright (C) 2009-2015  Syed Asad Rahman <asad @ ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.filters;

import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;

/**
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 * @author maclean
 * 
 */
public class BaseFilter {
    private static final ILoggingTool logger = createLoggingTool(BaseFilter.class);
    private static final Logger LOG = getLogger(BaseFilter.class.getName());

    private final IAtomContainer mol1;
    private final IAtomContainer mol2;

    /**
     *
     * @param sourceMol
     * @param targetMol
     */
    public BaseFilter(IAtomContainer sourceMol, IAtomContainer targetMol) {
        this.mol1 = sourceMol;
        this.mol2 = targetMol;

        try {
            percieveAtomTypesAndConfigureAtoms(mol1);
        } catch (CDKException ex) {
            logger.error(SEVERE, null, ex);
        }
        try {
            percieveAtomTypesAndConfigureAtoms(mol2);
        } catch (CDKException ex) {
            logger.error(SEVERE, null, ex);
        }
    }

    /**
     *
     * @param sourceMol
     * @param targetMol
     */
    public BaseFilter(IQueryAtomContainer sourceMol, IAtomContainer targetMol) {
        this.mol1 = sourceMol;
        this.mol2 = targetMol;

        try {
            percieveAtomTypesAndConfigureAtoms(mol2);
        } catch (CDKException ex) {
            logger.error(SEVERE, null, ex);
        }
    }

    /**
     * @return the mol1
     */
    public synchronized IAtomContainer getQuery() {
        return mol1;
    }

    /**
     * @return the mol2
     */
    public synchronized IAtomContainer getTarget() {
        return mol2;
    }
}
