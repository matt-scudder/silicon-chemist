/*
 * MX Cheminformatics Tools for Java
 *
 * Copyright (c) 2007-2009 Metamolecular, LLC
 *
 * http://metamolecular.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Copyright (C) 2009-2015  Syed Asad Rahman <asad @ ebi.ac.uk>
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
 *
 */
package org.openscience.smsd.ring;

import java.util.ArrayList;
import static java.util.Collections.reverse;
import static java.util.Collections.unmodifiableList;
import java.util.List;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IAtom;

/**
 *
 * 
 * 
 * @author Richard L. Apodaca <rapodaca at metamolecular.com> 2009-2015, 
 *         Syed Asad Rahman <asad @ ebi.ac.uk> 2009-2015
 */
public class PathEdge {
    private static final Logger LOG = getLogger(PathEdge.class.getName());

    private final List<IAtom> atoms;

    /**
     *
     * @param atoms
     */
    public PathEdge(List<IAtom> atoms) {
        this.atoms = atoms;
    }

    /**
     *
     * @return
     */
    public List<IAtom> getAtoms() {
        return unmodifiableList(atoms);
    }

    /**
     *
     * @return
     */
    public IAtom getSource() {
        return atoms.get(0);
    }

    /**
     *
     * @return
     */
    public IAtom getTarget() {
        return atoms.get(atoms.size() - 1);
    }

    /**
     *
     * @return
     */
    public boolean isCycle() {
        return (atoms.size() > 2) && atoms.get(0).equals(atoms.get(atoms.size() - 1));
    }

    /**
     *
     * @param other
     * @return
     */
    public PathEdge splice(PathEdge other) {
        IAtom intersection = getIntersection(other.atoms);
        List<IAtom> newAtoms = new ArrayList<>(atoms);

        if (atoms.get(0) == intersection) {
            reverse(newAtoms);
        }

        if (other.atoms.get(0) == intersection) {
            for (int i = 1; i < other.atoms.size(); i++) {
                newAtoms.add(other.atoms.get(i));
            }
        } else {
            for (int i = other.atoms.size() - 2; i >= 0; i--) {
                newAtoms.add(other.atoms.get(i));
            }
        }

        if (!isRealPath(newAtoms)) {
            return null;
        }

        return new PathEdge(newAtoms);
    }

    private boolean isRealPath(List<IAtom> atoms) {
        for (int i = 1; i < atoms.size() - 1; i++) {
            for (int j = 1; j < atoms.size() - 1; j++) {
                if (i == j) {
                    continue;
                }

                if (atoms.get(i) == atoms.get(j)) {
                    return false;
                }
            }
        }

        return true;
    }

    private IAtom getIntersection(List<IAtom> others) {
        if (atoms.get(atoms.size() - 1) == others.get(0) || atoms.get(atoms.size() - 1) == others.get(others.size() - 1)) {
            return atoms.get(atoms.size() - 1);
        }

        if (atoms.get(0) == others.get(0) || atoms.get(0) == others.get(others.size() - 1)) {
            return atoms.get(0);
        }

        throw new RuntimeException("Couldn't splice - no intersection.");
    }
}
