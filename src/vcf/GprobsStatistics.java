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

import blbutil.Const;
import java.text.DecimalFormat;
import java.util.Arrays;
import main.AlleleProbs;
import main.GenotypeValues;

/**
 * Class {@code GprobsStatistics} has methods for computing statistics
 * from posterior genotype probabilities.
 *
 * The squared correlation statistics computed by this class can be derived
 * using the methods found in Appendix 1 of
 * "Browning BL and Browning SR, Am J Hum Genet 2009;84(2):210-23".
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class GprobsStatistics {

    private final Marker marker;
    private final int nSamples;
    private final float[] alleleFreq;

    private float sumCall = 0;
    private float sumSquareCall = 0;
    private float sumExpected = 0;
    private float sumExpectedSquare = 0;
    private float sumSquareExpected= 0;
    private float sumCallExpected = 0;

    /**
     * Constructs a new {@code GprobsStatistics} instance from the
     * specified scaled genotype probabilities.
     * @param gv scaled sample posterior genotype probabilities
     * @param marker a marker index
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= gv.nMarkers()}
     * @throws NullPointerException if {@code gv == null}
     */
    public GprobsStatistics(GenotypeValues gv, int marker) {
        int nAlleles = gv.marker(marker).nAlleles();
        this.marker = gv.marker(marker);
        this.nSamples = gv.nSamples();
        this.alleleFreq = new float[nAlleles];
        float[] alProbs = new float[nAlleles];
        float[] gtProbs = new float[3];
        for (int j=0; j<this.nSamples; ++j) {
            setProbs(gv, marker, j, gtProbs, alProbs);
            for (int a=0; a<nAlleles; ++a) {
                alleleFreq[a] += alProbs[a];
            }
            int call = maxIndex(gtProbs);
            float exp = (gtProbs[1] + 2*gtProbs[2]);
            float expSquare = (gtProbs[1] + 4*gtProbs[2]);
            sumCall += call;
            sumSquareCall += call*call;
            sumExpected += exp;
            sumExpectedSquare += expSquare;
            sumSquareExpected += (exp*exp);
            sumCallExpected += (call*exp);
        }
        float sum = sum(alleleFreq);
        divideBy(alleleFreq, sum);
    }

    private static void setProbs(GenotypeValues gv, int marker, int sample,
            float[] gtProbs, float[] alProbs) {
        Arrays.fill(gtProbs, 0.0f);
        Arrays.fill(alProbs, 0.0f);
        int gt = 0;
        for (int a2=0; a2<alProbs.length; ++a2) {
            for (int a1=0; a1<=a2; ++a1) {
                float gprob = gv.value(marker, sample, gt++);
                alProbs[a1] += gprob;
                alProbs[a2] += gprob;
                if (a2==0) {
                    gtProbs[0] += gprob;
                }
                else if (a1==0) {
                    gtProbs[1] += gprob;
                }
                else {
                    gtProbs[2] += gprob;
                }
            }
        }
        float sum = sum(gtProbs);
        divideBy(gtProbs, sum);
        divideBy(alProbs, 2*sum);
    }

    /**
     * Constructs a new {@code GprobsStatistics} instance from the
     * specified allele probabilities.
     * @param alleleProbs the allele probabilities
     * @param marker a marker index
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= alProbs.nMarkers()}
     * @throws NullPointerException if {@code alProbs == null}
     */
    public GprobsStatistics(AlleleProbs alleleProbs, int marker) {
        int nAlleles = alleleProbs.marker(marker).nAlleles();
        this.marker = alleleProbs.marker(marker);
        this.nSamples = alleleProbs.nSamples();
        this.alleleFreq = new float[nAlleles];
        float[] alProbs = new float[nAlleles];
        float[] gtProbs = new float[3];
        for (int j=0; j<this.nSamples; ++j) {
            setProbs(alleleProbs, marker, j, gtProbs, alProbs);
            for (int a=0; a<nAlleles; ++a) {
                alleleFreq[a] += alProbs[a];
            }
            int call = maxIndex(gtProbs);
            float exp = (gtProbs[1] + 2*gtProbs[2]);
            float expSquare = (gtProbs[1] + 4*gtProbs[2]);
            sumCall += call;
            sumSquareCall += call*call;
            sumExpected += exp;
            sumExpectedSquare += expSquare;
            sumSquareExpected += (exp*exp);
            sumCallExpected += (call*exp);
        }
        float sum = sum(alleleFreq);
        divideBy(alleleFreq, sum);
    }

    private static void setProbs(AlleleProbs ap, int marker, int sample,
            float[] gtProbs, float[] alProbs) {
        Arrays.fill(gtProbs, 0.0f);
        Arrays.fill(alProbs, 0.0f);
        for (int a2=0; a2<alProbs.length; ++a2) {
            for (int a1=0; a1<=a2; ++a1) {
                float gprob = ap.gtProb(marker, sample, a1, a2);
                if (a1 != a2) {
                    gprob += ap.gtProb(marker, sample, a2, a1);
                }
                alProbs[a1] += gprob;
                alProbs[a2] += gprob;
                if (a2==0) {
                    gtProbs[0] += gprob;
                }
                else if (a1==0) {
                    gtProbs[1] += gprob;
                }
                else {
                    gtProbs[2] += gprob;
                }
            }
        }
        float sum = sum(gtProbs);
        divideBy(gtProbs, sum);
        divideBy(alProbs, 2*sum);
    }

    private static int maxIndex(float[] fa) {
        int maxIndex = 0;
        for (int j=1; j<fa.length; ++j) {
            if (fa[j]>fa[maxIndex]) {
                maxIndex = j;
            }
        }
        return maxIndex;
    }

    private static float sum(float[] fa) {
        float sum = 0.0f;
        for (float f : fa) {
            sum += f;
        }
        return sum;
    }

    private static void divideBy(float[] fa, float divisor) {
        for (int j=0; j<fa.length; ++j) {
            fa[j] /= divisor;
        }
    }

    /**
     * Returns the marker.
     * @return the marker.
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns an array of length {@code this.marker().nAlleles()} whose
     * {@code j}-th element is the estimated sample frequency of allele
     * {@code j}.
     * @return an array of length {@code this.marker().nAlleles()} whose
     * {@code j}-th element is the estimated sample frequency of allele
     * {@code j}
     */
    public float[] alleleFreq() {
        return alleleFreq.clone();
    }

    /**
     * Returns the estimated squared correlation between the most probable
     * ALT allele dose and the true ALT allele dose.
     * Returns 0 if the marker is monomorphic or if most probable ALT
     * allele dose is monomorphic.
     *
     * @return the estimated squared correlation between the most likely
     * allele dose and the true allele dose
     */
    public float allelicR2() {
        float f = 1.0f /  nSamples;
        float cov = sumCallExpected - (sumCall * sumExpected * f);
        float varBest = sumSquareCall - (sumCall * sumCall * f);
        float varExp = sumExpectedSquare - (sumExpected * sumExpected * f);
        float den = varBest * varExp;
        return (den==0.0f) ? 0.0f : Math.abs( (cov*cov) / den );
    }

    /**
     * Returns the estimated squared correlation between the estimated
     * ALT allele dose and the true ALT allele dose.  Returns 0 if the
     * marker is monomorphic.
     *
     * @return the estimated squared correlation between the estimated
     * ALT allele dose and the true ALT allele dose
     */
    public float doseR2() {
        float f = 1.0f / (float) nSamples;
        float num = sumSquareExpected - (sumExpected * sumExpected * f);
        float den = sumExpectedSquare - (sumExpected * sumExpected * f);
        return (den==0.0f) ? 0.0f : Math.abs(num / den);
    }

    /**
     * Returns the estimated squared correlation between the estimated
     * ALT allele dose and the true ALT allele dose where the variance of
     * the true ALT allele dose is estimated from the estimated
     * ALT allele frequency. Returns 0 if the marker is monomorphic.
     *
     * @return the estimated squared correlation between the estimated
     * ALT allele dose and the true ALT allele dose
     */
    public float hweDoseR2() {
        float f = 1.0f / nSamples;
        float num = (sumSquareExpected - (sumExpected*sumExpected*f))/nSamples;
        float altFreq = sumExpected / (2.0f * nSamples);
        float den = 2.0f * altFreq * (1.0f - altFreq);
        return (den==0.0f) ? 0.0f : Math.abs(num / den);
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        DecimalFormat df = new DecimalFormat("0.####");
        StringBuilder sb = new StringBuilder(80);
        sb.append(marker);
        sb.append(Const.tab);
        for (int j=0; j<alleleFreq.length; ++j) {
            sb.append( (j==0) ? "AF=" : Const.comma);
            sb.append(alleleFreq[j]);
        }
        sb.append(Const.tab);
        sb.append("AR2=");
        sb.append(format(df, allelicR2()));
        sb.append(Const.tab);
        sb.append("DR2=");
        sb.append(format(df, doseR2()));
        sb.append(Const.tab);
        sb.append("HDR2=");
        sb.append(format(df, hweDoseR2()));
        return sb.toString();
    }

    private static String format(DecimalFormat df, float d) {
        if (Double.isNaN(d)) {
            return "NaN";
        }
        else if (d==Double.POSITIVE_INFINITY) {
            return "Infinity";
        }
        else if (d==Double.NEGATIVE_INFINITY) {
            return "-Infinity";
        }
        else {
            return df.format(d);
        }
    }
}
