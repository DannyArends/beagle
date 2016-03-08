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
package main;

import beagleutil.Samples;
import blbutil.Const;
import blbutil.FileUtil;
import blbutil.IntPair;
import ibd.IbdSegment;
import java.io.Closeable;
import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import vcf.VcfWriter;

/**
 * <p>Class {@code WindowWriter} writes VCF and IBD output data.
 * </p>
 * <p>Instances of class {@code WindowWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class WindowWriter implements Closeable {

    private static final DecimalFormat df2 = new DecimalFormat("#.##");

    private boolean isClosed = false;
    private boolean appendIbd = false;

    private final Samples samples;
    private final File vcfOutFile;
    private final File ibdOutFile;
    private final File hbdOutFile;
    private final PrintWriter vcfOut;
    private final Map<IntPair, IbdSegment> ibdBuffer = new HashMap<>();

    /**
     * Constructs a new {@code WindowWriter} object.
     * @param samples the sample whose data will be printed
     * @param outPrefix the output file prefix
     *
     * @throws IllegalArgumentException if {@code outPrefix.length() == 0}
     * @throws NullPointerException if
     * {@code samples == null || outPrefix == null}
     */
    public WindowWriter(Samples samples, String outPrefix) {
        if (samples==null) {
            throw new NullPointerException("samples==null");
        }
        if (outPrefix.length()==0) {
            throw new IllegalArgumentException("outPrefix.length()==0");
        }
        this.samples = samples;
        this.vcfOutFile = new File(outPrefix + ".vcf.gz");
        this.ibdOutFile = new File(outPrefix + ".ibd");
        this.hbdOutFile = new File(outPrefix + ".hbd");
        this.vcfOut = FileUtil.bgzipPrintWriter(vcfOutFile);

        boolean printGT = true;
        boolean printGP = true;
        boolean printGL = false;
        VcfWriter.writeMetaLines(samples.ids(), Main.program,
                printGT, printGP, printGL, vcfOut);
    }

    /**
     * Returns the samples whose data is written by {@code this}.
     * @return the samples whose data is written by {@code this}
     */
    public Samples samples() {
        return samples;
    }


    /**
     * Returns {@code true} if {@code this.close()} method has
     * been previously invoked and returns {@code false} otherwise.
     *
     * @return {@code true} if {@code this.close()} method has
     * been previously invoked
     */
    public boolean isClosed() {
        return isClosed;
    }

    /**
     * Closes this {@code WindowWriter} for writing.  Calling the
     * {@code print()} method after invoking {@code close()} will
     * throw an {@code IllegalStateException}.
     */
    @Override
    public void close() {
        vcfOut.close();
        isClosed = true;
    }

    /**
     * Prints VCF records with GT and GP format fields for markers with
     * index between {@code cd.lastSplice()} (inclusive) and
     * {@code cd.nextSplice()} (exclusive).
     *
     * @param cd the input data for the current marker window
     * @param gv scaled genotype probabilities for the target samples
     *
     * @throws NullPointerException if {@code cd == null || gv == null}
     */
    public void printGV(CurrentData cd, GenotypeValues gv) {
        if (isClosed) {
            throw new IllegalStateException("isClosed()==true");
        }
        VcfWriter.appendRecords(gv, cd.prevTargetSplice(),
                cd.nextTargetSplice(), vcfOut);
        vcfOut.flush();
    }

    /**
     * Prints the data in {@code alProbs} for markers
     * with index between {@code cd.lastSplice()} (inclusive) and
     * {@code cd.nextSplice()} (exclusive).
     *
     * @param cd the input data for the current marker window
     * @param alProbs the estimated haplotype allele probabilities
     * @param imputed {@code true} if there are imputed markers,
     * and {@code false} otherwise
     * @param gprobs {@code true} if the GP field should be printed, and
     * {@code false} otherwise
     *
     * @throws IllegalStateException if {@code this.isClosed() == true}
     * @throws IllegalArgumentException if
     * {@code this.samples().equals(cd.targetSamples()) == false}
     * @throws IllegalArgumentException if
     * {@code this.samples().equals(alProbs.samples()) == false}
     * @throws IllegalArgumentException if
     * {@code cd.markers().equals(alProbs.markers()) == false}
     * @throws NullPointerException if {@code cd == null || alProbs == null}
     */
    public void print(CurrentData cd, AlleleProbs alProbs, boolean imputed,
            boolean gprobs) {
        if (isClosed) {
            throw new IllegalStateException("isClosed()==true");
        }
        if (cd.markers().equals(alProbs.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (samples.equals(cd.targetSamples()) == false
                || samples.equals(alProbs.samples()) == false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        int start = cd.prevSplice();
        int end = cd.nextSplice();
        VcfWriter.appendRecords(alProbs, start, end, imputed, gprobs, vcfOut);
        vcfOut.flush();
    }

    /**
     * Prints IBD segments that end between the markers
     * with index between {@code cd.lastSplice()} (inclusive) and
     * {@code cd.nextSplice()} (exclusive).
     * IBD segments that end on or after the marker with index
     * {@code cd.nextSplice()} are saved so that they can be merged
     * with IBD segments from the next marker window.
     *
     * <p>It is the the caller's responsibility to ensure that the ordered
     * haplotype pairs between adjacent consecutive markers windows
     * are identical for each sample.
     * </p>
     *
     * @param cd the input data for the current window
     * @param ibdMap a map whose keys are pairs of haplotype indices and whose
     * values are lists of IBD segments involving the haplotype pair key
     *
     * @throws IllegalStateException if {@code this.isClosed()==true}
     * @throws IllegalArgumentException if
     * {@code this.samples().equals(cd.targetSamples()) == false}
     * @throws NullPointerException if {@code cd == null || ibdMap == null}
     */
    public void printIbd(CurrentData cd, Map<IntPair, List<IbdSegment>> ibdMap) {
        if (isClosed) {
            throw new IllegalStateException("isClosed()==true");
        }
        if (samples.equals(cd.targetSamples()) == false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        printIbd(ibdMap, cd.prevTargetSplice(), cd.nextTargetOverlap(),
                cd.nextTargetSplice(), cd.nTargetMarkers());
        if (appendIbd==false) {
            appendIbd = true;
        }
    }

    private void printIbd(Map<IntPair, List<IbdSegment>> ibd, int lastSplice,
            int nextOverlap, int nextSplice, int nMarkers) {
        Map<IntPair, IbdSegment> lastBuffer = new HashMap<>(ibdBuffer);
        ibdBuffer.clear();
        try (PrintWriter ibdOut = FileUtil.printWriter(ibdOutFile, appendIbd);
               PrintWriter hbdOut = FileUtil.printWriter(hbdOutFile, appendIbd)) {
            Iterator<IntPair> keyIt = ibd.keySet().iterator();
            while (keyIt.hasNext()) {
                IntPair key = keyIt.next();
                List<IbdSegment> list = ibd.get(key);
                for (IbdSegment seg : list) {
                    if (seg.startIndex()==0) {
                        IbdSegment saved = lastBuffer.get(key);
                        if (saved!=null) {
                            seg = merge(saved, seg);
                        }
                    }
                    int ep1 = seg.endIndex()+1;
                    if (ep1>=lastSplice && (nextSplice==nMarkers || ep1<nextSplice)) {
                        printSegment(samples, seg, ibdOut, hbdOut);
                    }
                    else if (seg.startIndex()<nextOverlap) {
                        ibdBuffer.put(key, seg);
                    }
                }
                keyIt.remove();
            }
        }
    }

    private static IbdSegment merge(IbdSegment a, IbdSegment b) {
        assert a.hapPair().equals(b.hapPair());
        assert a.start().chromIndex()==b.start().chromIndex();
        int newStartIndex = -1;
        float newScore = Math.max(a.score(), b.score());
        return new IbdSegment(a.hapPair(), a.start(), b.end(),
                newScore, newStartIndex, b.endIndex());
    }

    private static void printSegment(Samples samples, IbdSegment tract,
            PrintWriter ibdOut, PrintWriter hbdOut) {
        int h1 = tract.hap1();
        int h2 = tract.hap2();
        int s1 = h1/2;
        int s2 = h2/2;
        PrintWriter out = (s1==s2) ? hbdOut : ibdOut;
        out.print(samples.id(s1));
        out.print(Const.tab);
        out.print((h1 % 2) + 1);
        out.print(Const.tab);
        out.print(samples.id(s2));
        out.print(Const.tab);
        out.print((h2 % 2) + 1);
        out.print(Const.tab);
        out.print(tract.start().chrom());
        out.print(Const.tab);
        out.print(tract.start().pos());
        out.print(Const.tab);
        out.print(tract.end().pos());
        out.print(Const.tab);
        out.println(df2.format(tract.score()));
    }
}
