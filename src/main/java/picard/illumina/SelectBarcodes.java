package picard.illumina;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.OtherProgramGroup;
import picard.util.StringDistanceUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Program to choose barcodes from a list of candidates and a distance requirement
 */
@CommandLineProgramProperties(
        summary = "...",
        oneLineSummary = "Program to choose barcodes from a list of candidates and a distance requirement",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
public class SelectBarcodes extends CommandLineProgram {

    @Argument
    public List<File> BARCODES_MUST_HAVE;

    @Argument
    public List<File> BARCODES_CHOOSE_FROM;

    @Argument
    public File SEED_BARCODES;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument
    public int FAR_ENOUGH = 3;

    @Argument
    public boolean COMPARE_REVCOMP = false;

    @Argument
    public boolean ALLOW_REVCOMP = false;

    @Argument
    public boolean ALLOW_REV = false;

    @Argument
    public boolean ALLOW_COMP = false;

    @Argument(optional = true)
    public File DISTANCES = null;

    private static final List<String> mustHaveBarcodes = new ArrayList<>();
    private static final List<String> barcodes = new ArrayList<>();
    private final List<String> seedBarcodes = new ArrayList<>();

    private final Map<Integer, BitSet> adjacencyMatrix = new HashMap<>();
    private final Map<Integer, BitSet> seedAdjacencyMatrix = new HashMap<>();

    static File output = null;
    private static Log LOG = Log.getInstance(SelectBarcodes.class);

    // The P verticies are the "Potential" verticies, verticies that may be chosen
    private static final Map<Integer, BitSet> Ps = new CollectionUtil.DefaultingMap<>((i) -> new BitSet(), true);

    // The X verticies are the "eXcluded" verticies, those that should not be examined
    private static final Map<Integer, BitSet> Xs = new CollectionUtil.DefaultingMap<>((i) -> new BitSet(), true);

    // the R verticies are the "Required" verticies, those that must be part of the clique
    private static final Map<Integer, BitSet> Rs = new CollectionUtil.DefaultingMap<>((i) -> new BitSet(), true);

    // this is a temporary variable holding the difference between P and the neighbors of the "pivot" vertex
    private static final Map<Integer, BitSet> Diffs = new HashMap<>();
    private static int recursionLevel;



    @Override
    public int doWork() {
        output = OUTPUT;

        //open files make lists of barcodes
        openBarcodes();

        LOG.info("Opened Barcode files.");
        // calculate distance matrix and adjacency matrix

        calculateAdjacencyMatrix();
        LOG.info("Calculated distances");

        LOG.info("there are " + mustHaveBarcodes.size() + " MUST_HAVE barcodes.");
        LOG.info("there are " + barcodes.size() + " other barcodes to choose from (after possibly rejecting some).");

        try (final PrintWriter writer = new PrintWriter("all.barcodes.txt")) {
            mustHaveBarcodes.forEach(writer::println);
            barcodes.forEach(writer::println);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        final BitSet R = new BitSet();
        final BitSet P = new BitSet();
        final BitSet nodeSubset = new BitSet();

        // add each group nodes from the seeds and
        // call BronKerbosch2 saving best selections to disk.
        for (int i = 0; i < seedBarcodes.size(); i++) {
            LOG.info("Adding " + seedAdjacencyMatrix.get(i).cardinality() + " nodes from seed " + i);
            nodeSubset.or(seedAdjacencyMatrix.get(i));

            final Map<Integer, BitSet> subGraph = subsetGraph(adjacencyMatrix, nodeSubset);
            final BitSet seedSolution = find_cliques(subsetGraph(adjacencyMatrix, nodeSubset), R);
            R.or(seedSolution);
        }
        // finally, add all remaining nodes

        LOG.info("Adding " + (adjacencyMatrix.size()-nodeSubset.cardinality()) + " remaining nodes.");

        final BitSet solution = find_cliques(adjacencyMatrix, R);

        LOG.info("final solution has cardinality " + solution.cardinality());

        return 0;
    }

    private void calculateAdjacencyMatrix() {

        final List<String> filteredBarcodes = barcodes.stream().filter(b -> {
                    final Optional<String> firstClose = mustHaveBarcodes
                            .stream()
                            .filter(m -> !areFarEnough(b, m))
                            .findAny();
                    if (firstClose.isPresent()) {

                        LOG.info(String.format("rejecting barcode: %s, it's too close to a MUST_HAVE barcode: %s.",
                                b, firstClose.get()));
                        return false;
                    }
                    return true;
                }
        ).collect(Collectors.toList());

        barcodes.clear();
        barcodes.addAll(filteredBarcodes);

        for (int ii = 0; ii < barcodes.size(); ii++) {
            final BitSet adjacency = new BitSet(barcodes.size());

            for (int jj = 0; jj < barcodes.size(); jj++) {
                adjacency.set(jj, areFarEnough(barcodes.get(ii), barcodes.get(jj)));
            }

            adjacencyMatrix.put(ii, adjacency);
        }

        if (DISTANCES != null) {
            try (final PrintWriter writer = new PrintWriter(DISTANCES)) {
                writer.append("BARCODE\t");
                writer.println(String.join("\t", barcodes));
                for (int ii = 0; ii < barcodes.size(); ii++) {
                    final BitSet adjacency = adjacencyMatrix.get(ii);

                    writer.append(barcodes.get(ii)).append('\t');
                    for (int jj = 0; jj < barcodes.size(); jj++) {
                        writer.append(adjacency.get(jj) ? "1" : "0").append('\t');
                    }
                    writer.append('\n');
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }


        for (int ii = 0; ii < seedBarcodes.size(); ii++) {
            final BitSet adjacency = new BitSet(barcodes.size());

            for (int jj = 0; jj < barcodes.size(); jj++) {
                adjacency.set(jj, areFarEnough(seedBarcodes.get(ii), barcodes.get(jj)));
            }

            seedAdjacencyMatrix.put(ii, adjacency);
        }
    }

    static int levenshtein(final String lhs, final String rhs, final int farEnough) {
        return StringDistanceUtils.levenshteinDistance(lhs.getBytes(), rhs.getBytes(), farEnough);
    }

    private boolean areFarEnough(final String lhs, final String rhs) {
        if (!ALLOW_REV) {
            final byte[] rev = rhs.getBytes();
            SequenceUtil.reverse(rev, 0, rev.length);
            if (lhs.equals(Arrays.toString(rev))) {
                return false;
            }
        }

        if (!ALLOW_COMP) {
            final byte[] comp = SequenceUtil.reverseComplement(rhs).getBytes();
            SequenceUtil.reverse(comp, 0, comp.length);
            if (lhs.equals(Arrays.toString(comp))) {
                return false;
            }
        }

        if (!ALLOW_REVCOMP && lhs.equals(SequenceUtil.reverseComplement(rhs))) {
            return false;
        }

        return levenshtein(lhs, rhs, FAR_ENOUGH) >= FAR_ENOUGH &&
                (!COMPARE_REVCOMP || levenshtein(lhs, SequenceUtil.reverseComplement(rhs), FAR_ENOUGH) >= FAR_ENOUGH );
    }

    //returns the number of SEED barcodes that were placed into "barcodes"
    private void openBarcodes() {
        barcodes.clear();
        mustHaveBarcodes.clear();

        BARCODES_MUST_HAVE.forEach(b -> readBarcodesFile(b, mustHaveBarcodes));
        BARCODES_CHOOSE_FROM.forEach(b -> readBarcodesFile(b, barcodes));

        //shuffle input barcodes to prevent bias
        Collections.shuffle(barcodes, new Random(51));
    }

    private void readBarcodesFile(final File f, final List<String> addTo) {
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().equals("")) continue;
                if (!addTo.contains(line)) {
                    addTo.add(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /** implements BronKerbosch3 https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
     * to find the maximal clique in the graph defined by the edges in {@code graph}. The
     * original algorithm is defined recursively, but here the recursion has been unfolded explicitly
     * to avoid the recursion overhead and to reuse objects.
     *
     *      BronKerbosch3(G):
     *        P = V(G)
     *        R = X = empty
     *        for each vertex v in a degeneracy ordering of G:
     *            BronKerbosch2({v}, P ⋂ N(v), X ⋂ N(v))
     *            P := P \ {v}
     *            X := X ⋃ {v}
     *
     *
     *      BronKerbosch2(R,P,X):
     *        if P and X are both empty:
     *            report R as a maximal clique
     *        choose a pivot vertex u in P ⋃ X
     *        for each vertex v in P \ N(u):
     *            BronKerbosch2(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
     *            P := P \ {v}
     *            X := X ⋃ {v}
     *
     *
     * @param graph the input bit-sets defining the edges between barcodes that are compatible
     */
    static BitSet find_cliques(Map<Integer, BitSet> graph, BitSet required) {

        final Integer[] degeneracyOrder = getDegeneracyOrder(graph);
        recursionLevel = 0;


        final BitSet best_clique = new BitSet();
        int bestCliqueSize = 0;

        for (final Integer v : degeneracyOrder) {
            recursionLevel = 0;

            BitSet r, p, x;

            r = Rs.get(recursionLevel);
            r.clear();
            r.or(required);

            x = Xs.get(recursionLevel);
            x.clear();

            // any node that isn't a neighbor of all the required nodes should be removed from consideration
            graph.keySet().stream()
                    .filter(i -> intersection(graph.get(i), required).cardinality() != required.cardinality())
                    .forEach(x::set);

            if (r.get(v) || x.get(v)) continue;

            r.set(v);

            p = Ps.get(recursionLevel);
            LOG.info("examining node " + v);
            p.clear();
            p.or(graph.get(v));

            while (recursionLevel >= 0) {

                p = Ps.get(recursionLevel);
                r = Rs.get(recursionLevel);
                x = Xs.get(recursionLevel);

                //  LOG.info(String.format("in while (%d) %s %s %s (%s)", recursionLevel, r, p, x, best_clique));

                // if there are no more pivots to chose from, we have reached a (locally) maximal clique.
                // or if there is no way we could find a larger clique, stop trying
                if (p.isEmpty() || p.cardinality() + r.cardinality() <= bestCliqueSize) {

                    Diffs.remove(recursionLevel);
                    registerClique(r, best_clique);
                    bestCliqueSize = best_clique.cardinality();
                    recursionLevel--;

                    continue;
                }

                if (!Diffs.containsKey(recursionLevel)) {
                    final BitSet finalP = p;

                    final int u = Stream.concat(p.stream().boxed(), x.stream().boxed())
                            .max(Comparator.comparingInt(o -> intersection(finalP, graph.get(o)).cardinality()))
                            .orElse(-1); //should never happen

                    Diffs.put(recursionLevel, difference(p, graph.get(u)));
                }

                final int vv = Diffs.get(recursionLevel).nextSetBit(0);
                if (vv == -1) {
                    Diffs.remove(recursionLevel);
                    p.clear();
                    x.clear();
                    r.clear();

                    recursionLevel--;
                    continue;
                }

                Diffs.get(recursionLevel).clear(vv);

                final BitSet recNeighs = graph.get(vv);

                recursionLevel++;

                final BitSet lowerR = Rs.get(recursionLevel);
                lowerR.clear();
                lowerR.or(r);
                lowerR.set(vv);

                final BitSet lowerX = Xs.get(recursionLevel);
                lowerX.clear();
                lowerX.or(x);
                lowerX.and(recNeighs);

                final BitSet lowerP = Ps.get(recursionLevel);
                lowerP.clear();
                lowerP.or(p);
                lowerP.and(recNeighs);

                p.clear(vv);
                x.set(vv);

            }

            p.clear(v);
            x.set(v);
        }
        return best_clique;
    }

    private static void registerClique(final BitSet r, final BitSet bestClique) {
        if (r.cardinality() > bestClique.cardinality()) {
            bestClique.clear();
            bestClique.or(r);
            System.out.println("best.cardinality()=" + (bestClique.cardinality() + mustHaveBarcodes.size()));

            output.delete();
            try (PrintWriter writer = new PrintWriter(output)) {

                mustHaveBarcodes.forEach(b -> writer.println(-1 + "\t" + b));
                bestClique.stream().forEach(i -> writer.println(i + "\t" + barcodes.get(i)));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    static Map<Integer,BitSet> subsetGraph(final Map<Integer,BitSet> graph, final BitSet mask){

        final Map<Integer, BitSet> retVal = new HashMap<>();

        mask.stream().forEach(i-> {
            retVal.put(i, intersection(graph.get(i),mask));
        });

        return retVal;
    }

    static BitSet union(final BitSet lhs, final BitSet rhs) {
        BitSet ret = BitSet.valueOf(lhs.toLongArray());
        ret.or(rhs);
        return ret;
    }

    static BitSet intersection(final BitSet lhs, final BitSet rhs) {
        BitSet ret = BitSet.valueOf(lhs.toLongArray());
        ret.and(rhs);
        return ret;
    }

    static BitSet difference(final BitSet lhs, final BitSet rhs) {
        BitSet ret = BitSet.valueOf(lhs.toLongArray());
        ret.andNot(rhs);
        return ret;
    }

    private static Integer[] getDegeneracyOrder(final Map<Integer, BitSet> graph) {
        List<Integer> ordering = new ArrayList<>();
        Set<Integer> ordering_set = new HashSet<>();

        // a map from the vertices to their cardinality
        Map<Integer, Integer> degrees = new CollectionUtil.DefaultingMap<>(0);

        // a map form a given degeneracy to the list of vertices with that degeneracy
        @SuppressWarnings("MismatchedQueryAndUpdateOfCollection")
        Map<Integer, List<Integer>> degen = new CollectionUtil.DefaultingMap<>(i -> new ArrayList<>(), true);
        int max_deg = -1;
        for (int v = 0; v < graph.size(); v++) {
            int deg = graph.get(v).cardinality();
            degen.get(deg).add(v);
            degrees.put(v, deg);
            if (deg > max_deg) {
                max_deg = deg;
            }
        }
        outter:
        while (true) {
            int i = 0;
            while (true) {
                if (i <= max_deg) {
                    if (degen.get(i).size() != 0) {
                        break;
                    }
                    i += 1;

                } else {
                    break outter;
                }
            }

            final Integer v = degen.get(i).remove(degen.get(i).size() - 1);

            ordering.add(v);
            ordering_set.add(v);
            graph.get(v).stream().forEach(w -> {

                if (!ordering_set.contains(w)) {
                    final int deg = degrees.get(w);
                    degen.get(deg).remove(Integer.valueOf(w));
                    if (deg > 0) {
                        degrees.put(w, degrees.get(w) - 1);
                        degen.get(deg - 1).add(w);
                    }
                }
            });
        }
        Collections.reverse(ordering);
        return ordering.toArray(new Integer[0]);
    }
}