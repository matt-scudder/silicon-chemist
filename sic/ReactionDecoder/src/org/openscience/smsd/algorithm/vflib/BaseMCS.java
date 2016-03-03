/* 
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
 * You should have received commonAtomList copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.vflib;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import static java.util.Collections.sort;
import static java.util.Collections.synchronizedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.mcgregor.McGregor;
import static org.openscience.smsd.algorithm.vflib.SortOrder.DESCENDING;
import org.openscience.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.smsd.algorithm.vflib.interfaces.IQuery;

/**
 * This class should be used to find MCS between source graph and target graph.
 *
 * First the algorithm runs VF lib
 * {@link org.openscience.cdk.smsd.algorithm.vflib.map.VFMCSMapper} and reports
 * MCS between run source and target graphs. Then these solutions are extended
 * using McGregor {@link org.openscience.cdk.smsd.algorithm.mcgregor.McGregor}
 * algorithm where ever required.
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BaseMCS {
    private static final ILoggingTool Logger = createLoggingTool(BaseMCS.class);

    /**
     *
     */
    protected int countR;

    /**
     *
     */
    protected int countP;

    /**
     *
     */
    protected final IAtomContainer source;

    /**
     *
     */
    protected final IAtomContainer target;
    private final boolean shouldMatchRings;
    private final boolean matchBonds;
    private final boolean matchAtomType;

    /**
     *
     */
    protected final List<Map<INode, IAtom>> vfLibSolutions;
    final List<Map<Integer, Integer>> allLocalMCS;
    final List<AtomAtomMapping> allLocalAtomAtomMapping;

    BaseMCS(IAtomContainer source, IAtomContainer target, boolean matchBonds, boolean shouldMatchRings, boolean matchAtomType) {
        this.allLocalAtomAtomMapping = new ArrayList<>();
        this.allLocalMCS = new ArrayList<>();
        this.shouldMatchRings = shouldMatchRings;
        this.matchBonds = matchBonds;
        this.matchAtomType = matchAtomType;
        this.vfLibSolutions = new ArrayList<>();
        this.source = source;
        this.target = target;
    }

    BaseMCS(IQueryAtomContainer source, IAtomContainer target) {
        this.allLocalAtomAtomMapping = new ArrayList<>();
        this.allLocalMCS = new ArrayList<>();
        this.shouldMatchRings = true;
        this.matchBonds = true;
        this.matchAtomType = true;
        this.vfLibSolutions = new ArrayList<>();
        this.source = source;
        this.target = target;
    }

    /**
     *
     * @param cliqueMap
     * @param mapGlobal
     * @return
     */
    protected synchronized boolean hasClique(
            Map<Integer, Integer> cliqueMap, List<Map<Integer, Integer>> mapGlobal) {
        for (Map<Integer, Integer> storedMap : mapGlobal) {
            if (cliqueMap.size() < storedMap.size()) {
                return true;
            } else if (cliqueMap.equals(storedMap)) {
                return true;
            }
        }
        return false;
    }

    /**
     *
     * @param cliqueMap
     * @param mapGlobal
     * @return
     */
    protected synchronized boolean isCliquePresent(Map<Integer, Integer> cliqueMap, Collection<Map<Integer, Integer>> mapGlobal) {
        for (Map<Integer, Integer> storedMap : mapGlobal) {
            if (cliqueMap.equals(storedMap)) {
                return true;
            }
        }
        return false;
    }

    /**
     *
     * @param refinedMCSSeeds
     * @throws CDKException
     * @throws IOException
     */
    protected synchronized void extendCliquesWithMcGregor(
            List<Map<Integer, Integer>> refinedMCSSeeds) throws CDKException, IOException {
        List<List<Integer>> mappings = new ArrayList<>();
        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : refinedMCSSeeds) {
            Map<Integer, Integer> extendMapping = new TreeMap<>(firstPassMappings);
            McGregor mgit;
            if (source instanceof IQueryAtomContainer) {
                mgit = new McGregor(source, target, mappings, isBondMatchFlag(), isMatchRings(), isMatchAtomType());
                //Start McGregor search
                mgit.startMcGregorIteration(source, mgit.getMCSSize(), extendMapping);
            } else {
                if (countR > countP) {
                    mgit = new McGregor(source, target, mappings, isBondMatchFlag(), isMatchRings(), isMatchAtomType());

                    //Start McGregor search
                    mgit.startMcGregorIteration(source, mgit.getMCSSize(), extendMapping);
                } else {
                    extendMapping.clear();
                    mgit = new McGregor(target, source, mappings, isBondMatchFlag(), isMatchRings(), isMatchAtomType());
                    ROPFlag = false;
                    for (Map.Entry<Integer, Integer> map : firstPassMappings.entrySet()) {
                        extendMapping.put(map.getValue(), map.getKey());
                    }
                    //Start McGregor search
                    mgit.startMcGregorIteration(target, mgit.getMCSSize(), extendMapping);
                }
            }
            mappings = mgit.getMappings();
        }
//        System.out.println("\nSol count after MG " + mappings.size());
        setMcGregorMappings(ROPFlag, mappings);
//        System.out.println("After set Sol count MG" + allMCS.size());
//        System.out.println("MCSSize " + vfMCSSize + "\n");
    }

    /**
     *
     * @param RONP
     * @param query
     */
    protected synchronized void setVFMappings(boolean RONP, IQuery query) {
        /*
         * Sort biggest clique to smallest
         */
        sort(vfLibSolutions, new Map2ValueComparator(DESCENDING));
        for (Map<INode, IAtom> solution : vfLibSolutions) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<>();

            for (INode node : solution.keySet()) {
                IAtom qAtom;
                IAtom tAtom;
                int qIndex;
                int tIndex;

                if (RONP) {
                    qAtom = query.getAtom(node);
                    tAtom = solution.get(node);
                    qIndex = source.getAtomNumber(qAtom);
                    tIndex = target.getAtomNumber(tAtom);
                } else {
                    tAtom = query.getAtom(node);
                    qAtom = solution.get(node);
                    qIndex = source.getAtomNumber(qAtom);
                    tIndex = target.getAtomNumber(tAtom);
                }

                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(SEVERE, null, ex);
                    }
                }
            }

            if (!indexindexMapping.isEmpty()
                    && !hasClique(indexindexMapping, getLocalMCSSolution())) {
                getLocalAtomMCSSolution().add(atomatomMapping);
                getLocalMCSSolution().add(indexindexMapping);
            }
        }
    }

    private synchronized void setMcGregorMappings(boolean RONP,
            List<List<Integer>> mappings) throws CDKException {
        int counter = 0;
        int solSize = 0;
        getLocalAtomMCSSolution().clear();
        getLocalMCSSolution().clear();
        for (List<Integer> mapping : mappings) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);
            Map<Integer, Integer> indexindexMapping = new TreeMap<>();
            for (int index = 0; index < mapping.size(); index += 2) {
                IAtom qAtom;
                IAtom tAtom;
                int qIndex;
                int tIndex;

                if (RONP) {
                    qAtom = source.getAtom(mapping.get(index));
                    tAtom = target.getAtom(mapping.get(index + 1));
                    qIndex = mapping.get(index);
                    tIndex = mapping.get(index + 1);
                } else {
                    qAtom = source.getAtom(mapping.get(index + 1));
                    tAtom = target.getAtom(mapping.get(index));
                    qIndex = mapping.get(index + 1);
                    tIndex = mapping.get(index);
                }

                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    throw new CDKException("Atom index pointing to NULL");
                }
            }
            if (indexindexMapping.size() > solSize) {
                solSize = indexindexMapping.size();
                getLocalAtomMCSSolution().clear();
                getLocalMCSSolution().clear();
                counter = 0;
            }
            if (!indexindexMapping.isEmpty()
                    && !hasClique(indexindexMapping, getLocalMCSSolution())
                    && indexindexMapping.size() == solSize) {
                getLocalAtomMCSSolution().add(counter, atomatomMapping);
                getLocalMCSSolution().add(counter, indexindexMapping);
                counter++;
            }
        }

    }

    /**
     *
     * @return
     */
    protected synchronized IAtomContainer getReactantMol() {
        return source;
    }

    /**
     *
     * @return
     */
    protected synchronized IAtomContainer getProductMol() {
        return target;
    }

    /**
     * @return the shouldMatchRings
     */
    protected boolean isMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @return the shouldMatchBonds
     */
    protected synchronized boolean isBondMatchFlag() {
        return matchBonds;
    }

    /**
     * @return the allLocalMCS
     */
    private synchronized List<Map<Integer, Integer>> getLocalMCSSolution() {
        return synchronizedList(allLocalMCS);
    }

    /**
     * @return the allLocalAtomAtomMapping
     */
    private synchronized List<AtomAtomMapping> getLocalAtomMCSSolution() {
        return synchronizedList(allLocalAtomAtomMapping);
    }

    /**
     *
     * @param mcsSeeds
     * @return
     */
    protected synchronized boolean isExtensionRequired(List<Map<Integer, Integer>> mcsSeeds) {
        int maxSize = 0;
        for (Map<Integer, Integer> map : mcsSeeds) {
            if (map.size() > maxSize) {
                maxSize = map.size();
            }
        }
        return this.source.getAtomCount() > maxSize && this.target.getAtomCount() > maxSize;
    }

    /**
     *
     * @return
     */
    protected synchronized boolean isExtensionRequired() {
        int commonAtomCount = checkCommonAtomCount(getReactantMol(), getProductMol());
        int maxSize = 0;
        for (Map<Integer, Integer> map : allLocalMCS) {
            if (map.size() > maxSize) {
                maxSize = map.size();
            }
        }
        return commonAtomCount > maxSize;
    }

    private synchronized int checkCommonAtomCount(
            IAtomContainer reactantMolecule, IAtomContainer productMolecule) {
        ArrayList<String> atoms = new ArrayList<>();
        for (int i = 0; i < reactantMolecule.getAtomCount(); i++) {
            atoms.add(reactantMolecule.getAtom(i).getSymbol());
        }
        int common = 0;
        for (int i = 0; i < productMolecule.getAtomCount(); i++) {
            String symbol = productMolecule.getAtom(i).getSymbol();
            if (atoms.contains(symbol)) {
                atoms.remove(symbol);
                common++;
            }
        }
        return common;
    }

    /**
     * @return the matchAtomType
     */
    public boolean isMatchAtomType() {
        return matchAtomType;
    }
    private static final java.util.logging.Logger LOG = java.util.logging.Logger.getLogger(BaseMCS.class.getName());
}
