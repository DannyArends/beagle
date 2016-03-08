/*
 * Copyright (C) 2014 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package sample;

import blbutil.IntPair;
import haplotype.SampleHapPairs;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * <p>Class {@code RefHapSegs} represents reference haplotypes that span
 * segments determined by non-overlapping clusters of markers.
 * </p>
 * <p>Instances of class {@code RefHapSegs} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefHapSegs {

    private final int[] clusterStart;
    private final int[] clusterEnd;
    private final SampleHapPairs refHapPairs;
    private final RefHapSeg[] refHapSegs;

    /**
     * Constructs a new {@code RefHapSegs} instance from the specified data.
     * @param refHapPairs the reference haplotype pairs
     * @param clusterStart an array whose {@code j}-th element is the
     * starting reference marker index (inclusive) for the {@code j}-th
     * marker cluster
     * @param clusterEnd an array whose {@code j}-th element is the
     * ending reference marker index (exclusive) for the {@code j}-th
     * marker cluster
     * @param nThreads the number of threads to use during object construction
     * @throws IllegalArgumentException if
     * {@code clusterStart.length != clusterEnd.length}
     * @throws IllegalArgumentException if
     * {@code clusterStart.length > 0 && clusterStart[0] < 0}
     * @throws IllegalArgumentException if
     * {@code clusterEnd.length > 0 && clusterEnd[clusterEnd.length - 1] > nMarkers}
     * @throws IllegalArgumentException if
     * {@code clusterStart[j] >= clusterEnd[j]} for some {@code j} satisfying
     * {@code 0 <= j && j < clusterStart.length}
     * @throws IllegalArgumentException if
     * {@code clusterStart[j] < clusterEnd[j-1]} for some {@code j} satisfying
     * {@code 1 <= j && j < clusterStart.length}
     * @throws IllegalArgumentException if {@code nThreads < 0}
     * @throws NullPointerException if
     * {@code refHapPairs == null || clusterStart == null || clusterEnd == null}
     */
    public RefHapSegs(SampleHapPairs refHapPairs, int[] clusterStart,
            int[] clusterEnd, int nThreads) {
        if (nThreads <=0) {
            throw new IllegalArgumentException(String.valueOf(nThreads));
        }
        int nMarkers = refHapPairs.nMarkers();
        checkClusters(clusterStart, clusterEnd, nMarkers);
        this.clusterStart = clusterStart.clone();
        this.clusterEnd = clusterEnd.clone();
        this.refHapPairs = refHapPairs;
        this.refHapSegs = IntStream.rangeClosed(0, this.clusterStart.length)
                .parallel()
                .mapToObj(j -> intPair(j, this.clusterStart, this.clusterEnd,
                        nMarkers))
                .map(ip -> new RefHapSeg(refHapPairs, ip.first(), ip.second()))
                .toArray(RefHapSeg[]::new);
    }

    private void checkClusters(int[] starts, int[] ends, int nMarkers) {
        if (starts.length != ends.length) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (starts.length > 0 && starts[0] < 0) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (ends.length > 0 && ends[ends.length - 1] > nMarkers) {
            throw new IllegalArgumentException("inconsistent data");
        }
        for (int j=0; j<starts.length; ++j) {
            if (starts[j] >= ends[j]) {
                throw new IllegalArgumentException("inconsistent data");
            }
            if (j>0 && ends[j-1] > starts[j]) {
                throw new IllegalArgumentException("inconsistent data");
            }
        }
    }

    private static IntPair intPair(int index, int[] starts, int[] ends,
            int nMarkers) {
        int start = (index == 0) ? 0 :  starts[index - 1];
        int end = (index == ends.length) ? nMarkers : ends[index];
        return new IntPair(start, end);
    }

    /**
     * Returns the reference haplotype pairs.
     * @return the reference haplotype pairs
     */
    public SampleHapPairs refHapPairs() {
        return refHapPairs;
    }

    /**
     * Return the number of distinct reference allele sequences in the
     * specified chromosome segment.
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @return the number of distinct reference allele sequences in the
     * specified chromosome segment
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     */
    public int nSeq(int segment) {
        return refHapSegs[segment].nSeq();
    }

    /**
     * Return the number of markers in the specified chromosome segment.
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @return the number of markers in the specified chromosome segment
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     */
    public int nMarkers(int segment) {
        return refHapSegs[segment].end() - refHapSegs[segment].start();
    }

    /**
     * Return the index of the allele sequence in the specified chromosome
     * segment for the specified reference haplotype.
     *
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @param hap a haplotype index
     *
     * @return the index of the allele sequence in the specified chromosome
     * segment for the specified reference haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.refHapPairs().nHaps()}
     */
    public int seq(int segment, int hap) {
        return refHapSegs[segment].seq(hap);
    }

    /**
     * Return the specified reference haplotype allele.
     *
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @param marker index of a marker in the specified interval
     * @param seq index of a reference allele sequence in the specified
     * interval
     * @return the specified reference haplotype allele
     *
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers(interval)}
     * @throws IndexOutOfBoundsException if
     * {@code seq < 0 || seg >= this.nSeq(segment)}
     */
    public int allele(int segment, int marker, int seq) {
        return refHapSegs[segment].allele(marker, seq);
    }

    /**
     * Returns the index of the first marker (inclusive) in the specified
     * marker cluster.
     * @param cluster an index of a marker cluster
     * @return the index of the first marker (inclusive) in the specified
     * marker cluster
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters}
     */
    public int clusterStart(int cluster) {
        return clusterStart[cluster];
    }

    /**
     * Returns the index of the last marker (exclusive) in the specified
     * marker cluster.
     * @param cluster an index of a marker cluster
     * @return the index of the last marker (exclusive) in the specified
     * marker cluster
     * @throws IndexOutOfBoundsException if
     * {@code cluster < 0 || cluster >= this.nClusters()}
     */
    public int clusterEnd(int cluster) {
        return clusterEnd[cluster];
    }

    /**
     * Returns the number of marker clusters.
     * @return the number of marker clusters
     */
    public int nClusters() {
        return clusterStart.length;
    }
}
