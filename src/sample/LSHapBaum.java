/*
 * Copyright (C) 2015 Brian L. Browning
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

import java.util.Arrays;
import main.HapAlleleProbs;
import main.LowMemHapAlleleProbs;
import vcf.Markers;

/**
 * <p>Class {@code LSHapBaum} implements the Baum hidden Markov model
 * forward and backward algorithms for imputing missing alleles on a
 * target haplotype.
 * </p>
 * <p>Instances of class {@code LSHapBaum} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class LSHapBaum {

    private final ImputationData impData;
    private final boolean lowMem;
    private final int n;    // number of reference haplotypes
    private final Markers refMarkers;
    private final float[] alleleProbs;
    private final float[][] fwdVal;
    private final float[] bwdVal;
    private final float[] emBwdVal;
    private final int[] fwdValueIndex2Marker;

    private final RefHapSegs refHapSegs;
    private final float[][] fwdHapProbs;
    private final float[][] bwdHapProbs;

    private float emBwdValuesSum = 0f;

    private int windowIndex = -9999;
    private int arrayIndex = -9999;

    /**
     * Creates a {@code LSHapBaum} instance from the specified data.
     *
     * @param impData the input data for genotype imputation
     * @param lowMem {@code true} if a low-memory checkpoint algorithm
     * should be used, and {@code false} otherwise
     *
     * @throws NullPointerException if {@code impData == null}
     */
    public LSHapBaum(ImputationData impData, boolean lowMem) {
        this.impData = impData;
        this.lowMem = lowMem;
        this.n = impData.refHapPairs().nHaps();
        this.refMarkers = impData.refHapPairs().markers();
        this.alleleProbs = new float[refMarkers.sumAlleles()];

        int nClusters = impData.nClusters();
        int size = lowMem ? (int) Math.ceil(Math.sqrt(1 + 8*nClusters)/2.0) + 1
                : nClusters;
        this.fwdValueIndex2Marker = new int[size];
        this.fwdVal = new float[size][n];
        this.bwdVal = new float[n];
        this.emBwdVal = new float[n];

        this.refHapSegs = impData.refHapSegs();
        this.fwdHapProbs = new float[impData.nClusters()][];
        this.bwdHapProbs = new float[impData.nClusters()][];
        for (int j=0; j < nClusters; ++j) {
            this.fwdHapProbs[j] = new float[refHapSegs.nSeq(j+1)];
            this.bwdHapProbs[j] = new float[refHapSegs.nSeq(j)];
        }
    }

    /**
     * <p>Estimates and returns allele probabilities for the specified target
     * haplotype. Estimated allele probabilities are conditional on the hidden
     * Markov model (HMM) and the input data represented by
     * {@code this.imputationData()}.
     * </p>
     *
     * @param hap a target data haplotype index
     * @return allele probabilities for the specified target haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.imputationData().targetHapPairs().nHaps()}
     */
    public HapAlleleProbs randomHapSample(int hap) {
        int nMarkers = impData.nClusters();
        Arrays.fill(alleleProbs, 0f);
        setFwdValues(hap);
        setInitBwdValue(hap);
        setStateProbs(nMarkers-1, currentIndex());
        for (int m=nMarkers-2; m>=0; --m) {
            setBwdValue(m, hap);
            setStateProbs(m, previousIndex(hap));
        }
        setAlleleProbs(alleleProbs);
        return new LowMemHapAlleleProbs(refMarkers, impData.targetSamples(),
                hap, alleleProbs);
    }

    /**
     * Returns the input data for genotype imputation.
     * @return the input data for genotype imputation
     */
    public ImputationData imputationData() {
        return impData;
    }

    private void setFwdValues(int hap) {
        int nMarkers = impData.nClusters();
        windowIndex = 0;
        arrayIndex = -1;
        for (int m=0; m<nMarkers; ++m) {
            float sum = 0f;
            float probRec = impData.pRecomb(m);
            int prev = currentIndex();
            int next = nextIndex();
            fwdValueIndex2Marker[next] = m;
            int a = impData.targetAllele(m, hap);
            for (int h=0; h<n; ++h) {
                int refAllele = impData.refAllele(m, h);
                float em = (a == refAllele) ? impData.noErrProb(m) : impData.errProb(m);
                float x = m==0 ? 1 : (probRec/n + (1-probRec)*fwdVal[prev][h]);
                fwdVal[next][h] = em*x;
                sum += fwdVal[next][h];
            }
            scale(fwdVal[next], sum);
        }
    }

    private static float sum(float[] fa) {
        float sum = 0f;
        for (float f : fa) {
            sum += f;
        }
        return sum;
    }

    private static void scale(float[] fa, float divisor) {
        for (int j=0; j<fa.length; ++j) {
            fa[j] /= divisor;
        }
    }

    private void setInitBwdValue(int hap) {
        int m = impData.nClusters() - 1;
        float f = 1f/n;
        emBwdValuesSum = 0f;
        int a = impData.targetAllele(m, hap);
        for (int h=0; h<n; ++h) {
            int refAllele = impData.refAllele(m, h);
            float em = (a == refAllele) ? impData.noErrProb(m) : impData.errProb(m);
            bwdVal[h] = f;
            emBwdVal[h] = f*em;
            emBwdValuesSum += emBwdVal[h];
        }
    }

    private void setBwdValue(int m, int hap) {
        float bwdValuesSum = 0f;
        float probRec = impData.pRecomb(m + 1);
        float commonTerm = emBwdValuesSum*probRec/n;
        for (int h=0; h<n; ++h) {
            bwdVal[h] = commonTerm + (1-probRec)*emBwdVal[h];
            bwdValuesSum += bwdVal[h];
        }
        int a = impData.targetAllele(m, hap);
        emBwdValuesSum = 0f;
        for (int h=0; h<n; ++h) {
            bwdVal[h] /= bwdValuesSum; // normalize first
            int refAllele = impData.refAllele(m, h);
            float em = (a == refAllele) ? impData.noErrProb(m) : impData.errProb(m);
            emBwdVal[h] = em*bwdVal[h];
            emBwdValuesSum += emBwdVal[h];
        }
    }

    private void setStateProbs(int m, int fwdIndex) {
        Arrays.fill(fwdHapProbs[m], 0f);
        Arrays.fill(bwdHapProbs[m], 0f);
        for (int h=0; h<n; ++h) {
            float stateProbs = fwdVal[fwdIndex][h]*bwdVal[h];
            fwdHapProbs[m][refHapSegs.seq(m+1, h)] += stateProbs;
            bwdHapProbs[m][refHapSegs.seq(m, h)] += stateProbs;
        }
        scale(fwdHapProbs[m], sum(fwdHapProbs[m]));
        scale(bwdHapProbs[m], sum(bwdHapProbs[m]));
    }

    private static float threshold(int nSeq) {
        return 0.5f/nSeq;
    }

    private void setAlleleProbs(float[] alleleProbs) {
        int nClusters = refHapSegs.nClusters();
        setFirstAlleleProbs(alleleProbs);
        for (int cluster=1; cluster < nClusters; ++cluster) {
            setAlleleProbs(alleleProbs, cluster);
        }
        setLastAlleleProbs(alleleProbs);
    }

    private void setFirstAlleleProbs(float[] alleleProbs) {
        int segment = 0;
        int refMarker = refHapSegs.clusterStart(segment);
        int nSeq = refHapSegs.nSeq(segment);
        float threshold = threshold(nSeq);
        for (int h=0; h<nSeq; ++h) {
            if (bwdHapProbs[segment][h] >= threshold) {
                for (int m=0; m<refMarker; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(segment, m, h);
                    alleleProbs[start + allele] += bwdHapProbs[segment][h];
                }
            }
        }
    }

    private void setAlleleProbs(float[] alleleProbs, int cluster) {
        assert cluster > 0;
        int startRefMarker = refHapSegs.clusterStart(cluster-1);
        int midRefMarker = refHapSegs.clusterEnd(cluster - 1);
        int endRefMarker = refHapSegs.clusterStart(cluster);
        int nSeq = refHapSegs.nSeq(cluster);
        float threshold = threshold(nSeq);
        for (int seq=0; seq<nSeq; ++seq) {
            boolean useFwd = fwdHapProbs[cluster-1][seq] >= threshold;
            boolean useBwd = bwdHapProbs[cluster][seq] >= threshold;
            if (useFwd) {
                for (int m=startRefMarker; m<midRefMarker; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(cluster, m - startRefMarker, seq);
                    alleleProbs[start + allele] += fwdHapProbs[cluster-1][seq];
                }
            }
            if (useFwd || useBwd) {
                for (int m=midRefMarker; m<endRefMarker; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(cluster, m - startRefMarker, seq);
                    double wt = impData.weight(m);
                    alleleProbs[start + allele] += wt*fwdHapProbs[cluster-1][seq];
                    alleleProbs[start + allele] += (1-wt)*bwdHapProbs[cluster][seq];
                }
            }
        }
    }

    private void setLastAlleleProbs(float[] alleleProbs) {
        int segment = refHapSegs.nClusters();
        int cluster = segment - 1;
        int refMarkerStart = refHapSegs.clusterStart(cluster);
        int refMarkerEnd = refHapSegs.refHapPairs().nMarkers();
        int nSeq = refHapSegs.nSeq(segment);
        float threshold = threshold(nSeq);
        for (int seq=0; seq<nSeq; ++seq) {
            if (fwdHapProbs[cluster][seq] >= threshold) {
                for (int m=refMarkerStart; m<refMarkerEnd; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(segment, m - refMarkerStart, seq);
                    alleleProbs[start + allele] += fwdHapProbs[cluster][seq];
                }
            }
        }
    }

    private int nextIndex() {
        ++arrayIndex;
        if (arrayIndex == fwdVal.length) {
            ++windowIndex;
            arrayIndex = windowIndex;
        }
        return arrayIndex;
    }

    private int currentIndex() {
        return arrayIndex;
    }

    private int previousIndex(int hap) {
        if (arrayIndex == windowIndex) {
            --windowIndex;
            arrayIndex = windowIndex;
            int start = fwdValueIndex2Marker[arrayIndex] + 1;
            int end = start + ( fwdVal.length - (arrayIndex + 1) );
            for (int m=start; m<end; ++m) {
                float sum = 0f;
                float probRec = impData.pRecomb(m);
                int prev = currentIndex();
                int next = nextIndex();
                fwdValueIndex2Marker[next] = m;
                int a = impData.targetAllele(m, hap);
                for (int h=0; h<n; ++h) {
                    int refAllele = impData.refAllele(m, h);
                    float em = (a == refAllele) ? impData.noErrProb(m) : impData.errProb(m);
                    float x = (probRec/n + (1-probRec)*fwdVal[prev][h]); // since m>0
                    fwdVal[next][h] = em*x;
                    sum += fwdVal[next][h];
                }
                scale(fwdVal[next], sum);
            }
            return arrayIndex;
        }
        else {
            return --arrayIndex;
        }
    }
}
