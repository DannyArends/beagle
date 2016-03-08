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
package vcf;

import beagleutil.Samples;
import blbutil.SampleFileIt;
import java.io.Closeable;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import main.GeneticMap;

/**
 * <p>Class {@code VcfWindow} represents a sliding window of VCF records.
 * </p>
 * Instances of class {@code VcfWindow} are not thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfWindow implements Closeable {

    private final SampleFileIt<? extends VcfEmission> it;
    private final List<VcfEmission> window;
    private int overlap;
    private int cumMarkerCnt;
    private VcfEmission next;

    /**
     * Constructs a new {@code VcfWindow} instance.
     * @param it an iterator that returns VCF records
     * @throws IllegalArgumentException if {@code it.hasNext() == false}
     * @throws IllegalArgumentException if a format error is detected in
     * a VCF record
     * @throws NullPointerException if {@code it == null}
     */
    public VcfWindow(SampleFileIt<? extends VcfEmission> it) {
        if (it.hasNext()==false) {
            throw new IllegalArgumentException("it.hasNext()==false");
        }
        this.it = it;
        this.overlap = 0;
        this.cumMarkerCnt = 0;
        this.window = new ArrayList<>(20000);
        this.next = it.next();
    }

    /**
     * Returns {@code true} if the sliding window of VCF Records is the last
     * window for the chromosome and returns {@code false} otherwise.
     * @return {@code true} if the sliding window of VCF Records is the last
     * window for the chromosome
     */
    public boolean lastWindowOnChrom() {
        return next==null || (sameChrom(next, window.get(0))==false);
    }

    private boolean sameChrom(VcfEmission a, VcfEmission b) {
        return a.marker().chromIndex()==b.marker().chromIndex();
    }

    /**
     * Returns {@code true} if the sliding window of VCF records can advance
     * and returns {@code false} otherwise.
     * @return {@code true} if the sliding window of VCF records can advance
     */
    public boolean canAdvanceWindow() {
        return next!=null;
    }

    /**
     * Advances the sliding window of VCF records, and returns the advanced
     * window as a {@code VcfEmission[]} object.  The size of the advanced
     * window and the number of markers of overlap between the marker window
     * immediately before method invocation and the marker window immediately
     * after method invocation may differ from the requested values.  If the
     * advanced window size or overlap is less than the requested value, the
     * actual value will be as large as possible. If
     * {@code this.lastWindowOnChrom() == true} before method invocation, then
     * there will be no overlap between the advanced window and the previous
     * window.
     *
     * @param overlap the requested number of markers of overlap
     * @param windowSize the requested number of the markers in the window
     * immediately after the method returns
     * @return the advanced window of VCF records
     *
     * @throws IllegalArgumentException if a format error is detected in
     * a VCF record
     * @throws IllegalArgumentException if
     * {@code overlap < 0 || overlap >= windowSize}
     * @throws IllegalStateException if
     * {@code this.canAdvanceWindow() == false}
     */
    public VcfEmission[] advanceWindow(int overlap, int windowSize) {
        if (canAdvanceWindow()==false) {
            throw new IllegalStateException("canAdvanceWindow()==false");
        }
        checkParameters(overlap, windowSize);
        overlap = getActualOverlap(overlap);
        List<VcfEmission> newWindow = new ArrayList<>(windowSize);

        newWindow.addAll(window.subList(window.size() - overlap, window.size()));
        int currentChromIndex = currentChromIndex(newWindow);
        while (newWindow.size() < windowSize
                && next != null
                && next.marker().chromIndex()==currentChromIndex) {
            newWindow.add(next);
            next = it.hasNext() ? it.next() : null;
        }
        // add all markers at the same marker position
        VcfEmission last = newWindow.get(newWindow.size()-1);
        while (next!=null && samePosition(last, next)) {
            newWindow.add(next);
            next = it.hasNext() ? it.next() : null;
        }
        this.overlap = overlap;
        this.window.clear();
        this.window.addAll(newWindow);
        this.cumMarkerCnt += (window.size() - overlap);
        return window.toArray(new VcfEmission[0]);
    }

    /**
     * Advances the sliding window of VCF records, and returns the advanced
     * window as a {@code VcfEmission[]} object.  The size of the advanced
     * window and the number of markers of overlap between the marker window
     * immediately before method invocation and the marker window immediately
     * after method invocation may differ from the requested values.  If the
     * distance the window is advanced or the overlap is less than the requested
     * value, the actual distance or overlap will be as large as possible. If
     * {@code this.lastWindowOnChrom() == true}
     * before method invocation, then there will be no overlap between the
     * advanced window and the previous window
     *
     * @param overlap the requested number of markers of overlap
     * @param cM the requested distance in cM to advance the window
     * @param map the genetic map
     * @return the advanced window of VCF records
     *
     * @throws IllegalArgumentException if a format error is detected in
     * a VCF record
     * @throws IllegalArgumentException if {@code overlap < 0 || cM <= 0}
     * @throws IllegalStateException if
     * {@code this.canAdvanceWindow() == false}
     */
    public VcfEmission[] advanceWindow(int overlap, double cM, GeneticMap map) {
        if (canAdvanceWindow()==false) {
            throw new IllegalStateException("canAdvanceWindow()==false");
        }
        if (overlap < 0) {
            throw new IllegalArgumentException(String.valueOf(overlap));
        }
        if (cM < 0) {
            throw new IllegalArgumentException(String.valueOf(cM));
        }

        overlap = getActualOverlap(overlap);
        List<VcfEmission> newWindow = new ArrayList<>(overlap + 1000);

        newWindow.addAll(window.subList(window.size() - overlap, window.size()));
        int currentChromIndex = currentChromIndex(newWindow);
        double endMapPos = startMapPos(newWindow, map) + cM;
        while (next != null
                && next.marker().chromIndex()==currentChromIndex
                && map.genPos(next.marker()) < endMapPos) {
            newWindow.add(next);
            next = it.hasNext() ? it.next() : null;
        }
        // add all markers at the same marker position
        VcfEmission last = newWindow.get(newWindow.size()-1);
        while (next!=null && samePosition(last, next)) {
            newWindow.add(next);
            next = it.hasNext() ? it.next() : null;
        }
        this.overlap = overlap;
        this.window.clear();
        this.window.addAll(newWindow);
        this.cumMarkerCnt += (window.size() - overlap);
        return window.toArray(new VcfEmission[0]);
    }

    private void checkParameters(int overlap, int windowSize) {
        if (overlap < 0 || overlap >= windowSize) {
            String s = "overlap=" + overlap + "windowSize=" + windowSize;
            throw new IllegalArgumentException(s);
        }
    }

    private int getActualOverlap(int overlap) {
        if (window.isEmpty() || lastWindowOnChrom()) {
            return 0;
        }
        int n = window.size();
        if (overlap > n) {
            overlap = n;
        }
        while (overlap > 0 && overlap < n
                && window.get(n - overlap).marker().pos()
                    == window.get(n - overlap - 1).marker().pos()) {
            ++overlap;
        }
        return overlap;
    }

    private int currentChromIndex(List<VcfEmission> currentWindow) {
        if (currentWindow.isEmpty()==false) {
            return currentWindow.get(0).marker().chromIndex();
        }
        else if (next!=null) {
            return next.marker().chromIndex();
        }
        else {
            return -1;
        }
    }

    private double startMapPos(List<VcfEmission> currentWindow, GeneticMap map) {
        if (currentWindow.isEmpty()==false) {
            Marker m = currentWindow.get(currentWindow.size() - 1).marker();
            return map.genPos(m);
        }
        else if (next!=null) {
            return map.genPos(next.marker());
        }
        else {
            return 0;
        }
    }

    private boolean samePosition(VcfEmission a, VcfEmission b) {
        return a.marker().chromIndex()==b.marker().chromIndex()
                && a.marker().pos()==b.marker().pos();
    }

    /**
     * Returns the file from which VCF records are read, or returns
     * {@code null} if the source is standard input.
     * @return the file from which VCF records are read, or
     * {@code null} if the source is standard input
     */
    public File file() {
        return it.file();
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return it.samples();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return  it.samples().nSamples();
    }

    /**
     * Returns the number of VCF records in the overlap between the current
     * window and the previous window.  Returns 0 if the current window
     * is the first window.
     *
     * @return the number of VCF records in the overlap between the current
     * window and the previous window
     */
    public int overlap() {
        return overlap;
    }

    /**
     * Returns the number of distinct VCF records in the union of the current
     * window and all previous windows.
     *
     * @return the number of distinct VCF records in the union of the current
     * window and all previous windows
     */
    public int cumMarkerCnt() {
        return cumMarkerCnt;
    }

    /**
     * Releases any I/O resources controlled by this object.
     */
    @Override
    public void close() {
        it.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(1100);
        sb.append(this.getClass().toString());
        sb.append("; next: ");
        sb.append(next);
        return sb.toString();
    }
}
