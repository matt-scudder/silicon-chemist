/*
 * Copyright (C) 2014 Syed Asad Rahman <asad at ebi.ac.uk>.
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
package org.openscience.smsd.algorithm.vflib.seeds;

import java.io.IOException;
import java.util.ArrayList;
import static java.util.Collections.sort;
import static java.util.Collections.unmodifiableList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.concurrent.Callable;
import static java.util.logging.Level.SEVERE;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.algorithm.mcsplus.BKKCKCF;
import org.openscience.smsd.algorithm.mcsplus.GenerateCompatibilityGraph;
import org.openscience.smsd.algorithm.rgraph.CDKRMapHandler;
import org.openscience.smsd.algorithm.vflib.Map1ValueComparator;
import static org.openscience.smsd.algorithm.vflib.SortOrder.DESCENDING;
import org.openscience.smsd.interfaces.Algorithm;
import static org.openscience.smsd.interfaces.Algorithm.CDKMCS;
import static org.openscience.smsd.interfaces.Algorithm.MCSPlus;

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
 * @author Syed Asad Rahman <asad at ebi.ac.uk>
 */
public class MCSSeedGenerator implements Callable<List<AtomAtomMapping>> {
    private static final ILoggingTool Logger = createLoggingTool(MCSSeedGenerator.class);

    private final IAtomContainer source;
    private final IAtomContainer target;
    private final List<AtomAtomMapping> allCliqueAtomMCS;
    private final boolean ringMatch;
    private final Algorithm algorithm;
    private final boolean bondMatch;
    private final boolean matchAtomType;

    /**
     *
     * @param source
     * @param target
     * @param bondMatch
     * @param ringMatch
     * @param matchAtomType
     * @param algorithm
     */
    public MCSSeedGenerator(IAtomContainer source, IAtomContainer target, boolean bondMatch, boolean ringMatch, boolean matchAtomType, Algorithm algorithm) {
        this.source = source;
        this.target = target;
        this.allCliqueAtomMCS = new ArrayList<>();
        this.ringMatch = ringMatch;
        this.algorithm = algorithm;
        this.matchAtomType = matchAtomType;
        this.bondMatch = bondMatch;
    }

    /**
     *
     * @param source
     * @param target
     * @param algorithm
     */
    public MCSSeedGenerator(IQueryAtomContainer source, IAtomContainer target, Algorithm algorithm) {
        this.source = source;
        this.target = target;
        this.allCliqueAtomMCS = new ArrayList<>();
        this.ringMatch = true;
        this.algorithm = algorithm;
        this.matchAtomType = true;
        this.bondMatch = true;
    }

    @Override
    public List<AtomAtomMapping> call() throws Exception {
//        System.out.println("ac1: " + this.source.getAtomCount());
//        System.out.println("ac2: " + this.target.getAtomCount());
        if (algorithm.equals(CDKMCS)) {
//            System.out.println("Calling CDKMCS " + bondMatch + " " + ringMatch);
            List<AtomAtomMapping> addUIT = addUIT();
//            System.out.println("addUIT " + addUIT.iterator().next().getCount());
            return addUIT;
        } else if (algorithm.equals(MCSPlus)) {
//            System.out.println("Calling MCSPLUS " + bondMatch + " " + ringMatch + " " + matchAtomType);
            List<AtomAtomMapping> addKochCliques = addKochCliques();
//            System.out.println("MCSPLUS " + addKochCliques.iterator().next().getCount());
            return addKochCliques;
        } else {
            return unmodifiableList(allCliqueAtomMCS);
        }
    }

    /**
     *
     * @return
     * @throws IOException
     */
    protected synchronized List<AtomAtomMapping> addKochCliques() throws IOException {
        IAtomContainer ac1;
        IAtomContainer ac2;
        boolean flagExchange = false;

        if (source instanceof IQueryAtomContainer) {
            ac1 = source;
            ac2 = target;
        } else if (source.getAtomCount() < target.getAtomCount()) {
            ac1 = source;
            ac2 = target;
        } else {
            flagExchange = true;
            ac1 = target;
            ac2 = source;
        }

        GenerateCompatibilityGraph gcg
                = new GenerateCompatibilityGraph(ac1, ac2, bondMatch, ringMatch, matchAtomType);
        List<Integer> comp_graph_nodes = gcg.getCompGraphNodes();
        List<Integer> cEdges = gcg.getCEgdes();
        List<Integer> dEdges = gcg.getDEgdes();
        BKKCKCF init = new BKKCKCF(comp_graph_nodes, cEdges, dEdges);
        Stack<List<Integer>> maxCliqueSet = new Stack<>();
        maxCliqueSet.addAll(init.getMaxCliqueSet());
        sort(maxCliqueSet, new Comparator<List<Integer>>() {
            @Override
            public int compare(List<Integer> a1, List<Integer> a2) {
                return a2.size() - a1.size(); // assumes you want biggest to smallest
            }
        });
        while (!maxCliqueSet.empty()) {
            List<Integer> peek = maxCliqueSet.peek();
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);

            for (Integer value : peek) {
                int[] index = getIndex(value, comp_graph_nodes);
                Integer qIndex = index[0];
                Integer tIndex = index[1];
                if (qIndex != -1 && tIndex != -1) {
                    IAtom qAtom;
                    IAtom tAtom;
                    if (flagExchange) {
                        qAtom = source.getAtom(tIndex);
                        tAtom = target.getAtom(qIndex);
                    } else {
                        qAtom = source.getAtom(qIndex);
                        tAtom = target.getAtom(tIndex);
                    }
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(SEVERE, null, ex);
                    }
                }
            }

            if (!atomatomMapping.isEmpty()) {
                allCliqueAtomMCS.add(atomatomMapping);
            }
            maxCliqueSet.pop();
        }
        gcg.clear();
        return unmodifiableList(allCliqueAtomMCS);
    }

    /**
     *
     * @return
     */
    private List<AtomAtomMapping> addUIT() throws CDKException {
        CDKRMapHandler rmap = new CDKRMapHandler();
        List<Map<Integer, Integer>> solutions;

        boolean rOnPFlag;
        if (source instanceof IQueryAtomContainer) {
            rOnPFlag = false;
            solutions = rmap.calculateOverlapsAndReduce(target, (IQueryAtomContainer) source);
        } else if (source.getAtomCount() > target.getAtomCount()) {
            rOnPFlag = true;
            solutions = rmap.calculateOverlapsAndReduce(source, target, bondMatch, ringMatch, matchAtomType);
        } else {
            rOnPFlag = false;
            solutions = rmap.calculateOverlapsAndReduce(target, source, bondMatch, ringMatch, matchAtomType);
        }
        return setUITMappings(rOnPFlag, solutions);
    }

    private List<AtomAtomMapping> setUITMappings(boolean RONP, List<Map<Integer, Integer>> sol) {
        /*
         * Sort biggest clique to smallest
         */
        sort(sol, new Map1ValueComparator(DESCENDING));
        for (Map<Integer, Integer> solution : sol) {
            AtomAtomMapping atomatomMapping = new AtomAtomMapping(source, target);

            for (Integer qAtomIndex : solution.keySet()) {
                IAtom qAtom;
                IAtom tAtom;
                int qIndex;
                int tIndex;

                if (RONP) {
                    qAtom = source.getAtom(qAtomIndex);
                    tAtom = target.getAtom(solution.get(qAtomIndex));
                } else {
                    tAtom = target.getAtom(qAtomIndex);
                    qAtom = source.getAtom(solution.get(qAtomIndex));
                }

                qIndex = source.getAtomNumber(qAtom);
                tIndex = target.getAtomNumber(tAtom);
                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to -1");
                    } catch (CDKException ex) {
                        Logger.error(SEVERE, null, ex);
                    }
                }
            }

            if (!atomatomMapping.isEmpty()) {
                allCliqueAtomMCS.add(atomatomMapping);
            }
        }
        return unmodifiableList(allCliqueAtomMCS);
    }

    private int[] getIndex(int cliqueIndex, List<Integer> comp_graph_nodes) {
        int[] v = new int[2];
        v[0] = -1;
        v[1] = -1;
        for (int i = 0; i < comp_graph_nodes.size(); i += 3) {
            if (cliqueIndex == comp_graph_nodes.get(i + 2)) {
                v[0] = comp_graph_nodes.get(i);
                v[1] = comp_graph_nodes.get(i + 1);
            }
        }
        return v;
    }
    private static final java.util.logging.Logger LOG = java.util.logging.Logger.getLogger(MCSSeedGenerator.class.getName());
}
