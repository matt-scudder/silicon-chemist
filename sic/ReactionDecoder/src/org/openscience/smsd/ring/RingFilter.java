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
import java.util.Collection;
import static java.util.Collections.sort;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.smsd.algorithm.matchers.AtomMatcher;
import static org.openscience.smsd.ring.HanserRingFinder.findRings;

/**
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk> 2009-2015
 */
public class RingFilter {
    private static final Logger LOG = getLogger(RingFilter.class.getName());

    private final Comparator comparator;
    private final AtomMatcher filter;
    private IAtomContainer mol;

    /**
     *
     * @param filter
     */
    public RingFilter(AtomMatcher filter) {
        comparator = new RingSizeComparator();
        this.filter = filter;
    }

    /**
     *
     * @param atomLimit
     * @param rings
     * @param atoms
     */
    public void filterAtoms(int atomLimit, Collection<List<IAtom>> rings, Collection<IAtom> atoms) {
        for (List<IAtom> ring : rings) {
            if (atoms.size() >= atomLimit) {
                break;
            }

            if (ringMatches(ring)) {
                atoms.addAll(ring);
            }
        }
    }

    /**
     *
     * @param molecule
     * @param atoms
     */
    public void filterAtoms(IAtomContainer molecule, Collection<IAtom> atoms) {
        this.mol = molecule;
        List<List<IAtom>> rings = new ArrayList<>(findRings(molecule));
        sort(rings, comparator);

        for (List<IAtom> ring : rings) {
            if (atoms.size() == molecule.getAtomCount()) {
                break;
            }

            if (ringMatches(ring)) {
                atoms.addAll(ring);
            }
        }
    }

    private boolean ringMatches(List<IAtom> ring) {
        for (IAtom atom : ring) {
            if (!filter.matches(atom)) {
                return false;
            }
        }

        return true;
    }

    private class RingSizeComparator implements Comparator<List<List<IAtom>>> {

        @Override
        public int compare(List<List<IAtom>> o1, List<List<IAtom>> o2) {
            if (o1.size() > o2.size()) {
                return +1;
            }
            if (o1.size() < o2.size()) {
                return -1;
            }
            return 0;
        }
    }
}
