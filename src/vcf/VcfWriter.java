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
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import main.AlleleProbs;
import main.GenotypeValues;

/**
 * <p>Class {@code VcfWriter} contains static methods for writing data in
 * VCF 4.2 format.
 * </p>
 * <p>Instances of class {@code VcfWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfWriter {

    private static final String PASS = "PASS";
    private static final DecimalFormat df2 = new DecimalFormat("#.##");
    private static final DecimalFormat df3 = new DecimalFormat("#.###");
    private static final DecimalFormat df2_fixed = new DecimalFormat("0.00");
    private static final MathContext mathContext2 = new MathContext(2);

    private static final String fileformat = "##fileformat=VCFv4.2";

    private static final String afInfo = "##INFO=<ID=AF,Number=A,Type=Float,"
            + "Description=\"Estimated Allele Frequencies\">";
    private static final String ar2Info = "##INFO=<ID=AR2,Number=1,Type=Float,"
            + "Description=\"Allelic R-Squared: estimated correlation between "
            + "most probable ALT dose and true ALT dose\">";
    private static final String dr2Info = "##INFO=<ID=DR2,Number=A,Type=Float,"
            + "Description=\"Dosage R-Squared: estimated correlation between "
            + "estimated ALT dose [P(RA) + 2*P(AA)] and true ALT dose\">";

    private static final String gtFormat = "##FORMAT=<ID=GT,Number=1,Type=String,"
            + "Description=\"Genotype\">";
    private static final String dsFormat = "##FORMAT=<ID=DS,Number=1,Type=Float,"
            +"Description=\"estimated ALT dose [P(RA) + P(AA)]\">";
    private static final String glFormat = "##FORMAT=<ID=GL,Number=G,Type=Float,"
            + "Description=\"Log10-scaled Genotype Likelihood\">";
    private static final String gpFormat = "##FORMAT=<ID=GP,Number=G,Type=Float,"
            + "Description=\"Estimated Genotype Probability\">";

    private static final String shortChromPrefix= "#CHROM" + Const.tab + "POS"
            + Const.tab + "ID" + Const.tab + "REF" + Const.tab + "ALT"
            + Const.tab + "QUAL" + Const.tab + "FILTER" + Const.tab + "INFO";

    private static final String longChromPrefix =
            shortChromPrefix + Const.tab + "FORMAT";


    private VcfWriter() {
        // private constructor prevents instantiation
    }

    /**
     * Writes VCF meta-information lines and header line to the specified
     * {@code PrintWriter}. Only one FORMAT subfield, the GT subfield,
     * is described in the meta-information lines.
     * @param sampleIds the sample identifiers
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed
     * @param out the {@code PrintWriter} to which VCF meta-information
     * lines will be written
     * @throws NullPointerException if {@code out == null}
     * @throws NullPointerException if
     * {@code sampleIds == null}, or if {@code sampleIds[j] == null} for any
     * {@code j} satisfying {@code (0 <= j && j < <sampleIds.length)}
     */
    public static void writeMetaLinesGT(String[] sampleIds, String source,
            PrintWriter out) {
        boolean printGT = true;
        boolean printGP = false;
        boolean printGL = false;
        writeMetaLines(sampleIds, source, printGT, printGP, printGL, out);
    }

    /**
     * Writes VCF meta-information lines and header line to the specified
     * {@code PrintWriter}.
     * @param sampleIds the sample identifiers
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed
     * @param printGT {@code true} if the meta-information lines
     * will describe the GT FORMAT subfield and {@code false} otherwise
     * @param printGP {@code true} if the meta-information lines
     * will describe the GP FORMAT subfield and {@code false} otherwise
     * @param printGL {@code true} if the meta-information lines
     * will describe the GL FORMAT subfield and {@code false} otherwise
     * @param out the {@code PrintWriter} to which VCF meta-information lines
     * will be written.
     * @throws NullPointerException if {@code out == null}
     * @throws NullPointerException if
     * {@code sampleIds == null}, or if {@code sampleIds[j] == null} for any
     * {@code j} satisfying {@code (0 <= j && j < sampleIds.length)}
     */
    public static void writeMetaLines(String[] sampleIds, String source,
            boolean printGT, boolean printGP, boolean printGL, PrintWriter out) {
        out.print(fileformat);
        out.print(Const.nl);
        out.print("##filedate=");
        out.print(now());
        out.print(Const.nl);
        if (source != null) {
            out.print("##source=\"");
            out.print(source);
            out.println("\"");
        }
        if (printGP) {
            out.println(afInfo);
            out.println(ar2Info);
            out.println(dr2Info);
        }
        if (printGT) {
            out.println(gtFormat);
        }
        if (printGL) {
            out.println(glFormat);
        }
        if (printGP) {
            out.println(dsFormat);
            out.println(gpFormat);
        }
        out.print(longChromPrefix);
        for (String id : sampleIds) {
            if (id==null) {
                throw new NullPointerException("id==null");
            }
            out.print(Const.tab);
            out.print(id);
        }
        out.println();
    }

    private static String now() {
        String dateFormat = "yyyyMMdd";
        Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat(dateFormat);
        return sdf.format(cal.getTime());
    }

    /**
     * Writes the specified genotype data  as VCF records to the specified
     * {@code PrintWriter}.
     * @param gv the scaled sample posterior genotype probabilities
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param out the {@code PrintWriter} to which VCF records will
     * be written.
     *
     * @throws IllegalArgumentException if
     * {@code haps.markers().equals(gv.markers()) == false}
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > haps.nMarkers())}
     * @throws NullPointerException if
     * {@code (gv == null || out == null)}
     */
    public static void appendRecords(GenotypeValues gv, int start, int end,
            PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        for (int marker=start; marker<end; ++marker) {
            printFixedFields(gv, marker, out);
            for (int sample=0, n=gv.nSamples(); sample<n; ++sample) {
                print_GT_DS_GP(gv, marker, sample, out);
            }
            out.println();
        }
    }


    private static void print_GT_DS_GP(GenotypeValues gv, int marker, int sample,
            PrintWriter out) {
        int nAlleles = gv.marker(marker).nAlleles();
        int nGenotypes = gv.marker(marker).nGenotypes();
        float[] dose = new float[nAlleles];
        int bestA1 = -1;
        int bestA2 = -1;
        int gt = 0;
        float sum = 0f;
        float maxGP = 0f;
        for (int a2=0; a2<nAlleles; ++a2) {
            for (int a1=0; a1<=a2; ++a1) {
                float value = gv.value(marker, sample, gt++);
                if (value > maxGP) {
                    bestA1 = a1;
                    bestA2 = a2;
                    maxGP = value;
                }
                dose[a1] += value;
                dose[a2] += value;
                sum += value;
            }
        }
        out.print(Const.tab);
        out.print(bestA1 == -1 ? Const.MISSING_DATA_STRING : bestA1);
        out.print(Const.unphasedSep);
        out.print(bestA2 == -1 ? Const.MISSING_DATA_STRING : bestA2);
        for (int al = 1; al < dose.length; ++al) {
            out.print( (al==1) ? Const.colon : Const.comma);
            out.print(df2.format(dose[al]/sum));
        }
        for (gt=0; gt<nGenotypes; ++gt) {
            out.print(gt==0 ? Const.colon : Const.comma);
            double v = gv.value(marker, sample, gt)/sum;
            out.print(df2.format(v));
        }
    }

    /**
     * Writes the specified genotype data as VCF records to the specified
     * {@code PrintWriter}.
     * @param alProbs the sample haplotype pairs
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param imputed {@code true} if there are imputed markers,
     * and {@code false} otherwise
     * @param gprobs {@code true} if the GP field should be printed, and
     * {@code false} otherwise.
     * @param out the {@code PrintWriter} to which VCF records will
     * be written
     *
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > alProbs.nMarkers())}
     * @throws NullPointerException if {@code haps == null || out == null}
     */
    public static void appendRecords(AlleleProbs alProbs, int start, int end,
            boolean imputed, boolean gprobs, PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        for (int marker=start; marker<end; ++marker) {
            printFixedFields(alProbs, marker, imputed, gprobs, out);
            for (int sample=0, n=alProbs.nSamples(); sample<n; ++sample) {
                printGTandDose(alProbs, marker, sample, imputed, out);
                if (gprobs) {
                     printGP(alProbs, marker, sample, out);
                }
            }
            out.println();
        }
    }

    private static void printGTandDose(AlleleProbs alProbs, int marker, int
            sample, boolean imputed, PrintWriter out) {
        out.print(Const.tab);
        out.print(alProbs.allele1(marker, sample));
        out.append(Const.phasedSep);
        out.print(alProbs.allele2(marker, sample));
        if (imputed) {
            int nAlleles = alProbs.marker(marker).nAlleles();
            for (int j = 1; j < nAlleles; ++j) {
                float p1 = alProbs.alProb1(marker, sample, j);
                float p2 = alProbs.alProb2(marker, sample, j);
                out.print( (j==1) ? Const.colon : Const.comma );
                out.print(df2.format(p1 + p2));
            }
        }
    }

    private static void printGP(AlleleProbs alProbs, int marker, int sample,
            PrintWriter out) {
        int nAlleles = alProbs.marker(marker).nAlleles();
        for (int a2=0; a2<nAlleles; ++a2) {
            for (int a1=0; a1<=a2; ++a1) {
                out.print((a2 == 0 && a1 == 0) ? Const.colon : Const.comma);
                float gtProb = alProbs.gtProb(marker, sample, a1, a2);
                if (a1 != a2) {
                    gtProb += alProbs.gtProb(marker, sample, a2, a1);
                }
                out.print(df2.format(gtProb));
            }
        }
    }

    /**
     * Prints the first 9 VCF record fields for the specified marker to
     * the specified {@code PrintWriter}.  Only one VCF FORMAT subfield,
     * the GT subfield, is printed.
     *
     * @param marker a marker
     * @param out the {@code PrintWriter} to which the first 9 VCF record
     * fields will be written
     *
     * @throws NullPointerException if {@code marker == null || out == null}
     */
    public static void printFixedFieldsGT(Marker marker, PrintWriter out) {
        out.print(marker);
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // INFO
        out.print(Const.tab);
        out.print("GT");
    }

    private static void printFixedFields(GenotypeValues gv, int marker,
            PrintWriter out) {
        GprobsStatistics gpm = new GprobsStatistics(gv, marker);
        float[] alleleFreq = gpm.alleleFreq();
        out.print(gv.marker(marker));
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        out.print(Const.tab);
        out.print("AR2=");                  // INFO
        out.print(df2_fixed.format(gpm.allelicR2()));
        out.print(";DR2=");
        out.print(df2_fixed.format(gpm.doseR2()));
        for (int j=1; j<alleleFreq.length; ++j) {
            out.print( (j==1) ? ";AF=" : Const.comma);
            BigDecimal bd = new BigDecimal(alleleFreq[j]).round(mathContext2);
            out.print(bd.doubleValue());
        }
        out.print(Const.tab);
        out.print("GT:DS:GP");
    }

    private static void printFixedFields(AlleleProbs alProbs,
            int marker, boolean printR2, boolean gprobs, PrintWriter out) {
        GprobsStatistics gpm = new GprobsStatistics(alProbs, marker);
        float[] alleleFreq = gpm.alleleFreq();
        out.print(alProbs.marker(marker));
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        if (printR2) {
            out.print(Const.tab);
            out.print("AR2=");                  // INFO
            out.print(df2_fixed.format(gpm.allelicR2()));
            out.print(";DR2=");
            out.print(df2_fixed.format(gpm.doseR2()));
            for (int j=1; j<alleleFreq.length; ++j) {
                out.print( (j==1) ? ";AF=" : Const.comma);
                BigDecimal bd = new BigDecimal(alleleFreq[j]).round(mathContext2);
                out.print(bd.doubleValue());
            }
        }
        else {
            out.print(Const.tab);
            out.print(Const.MISSING_DATA_CHAR);
        }
        out.print(Const.tab);
        if (printR2) {
            out.print(gprobs ? "GT:DS:GP" : "GT:DS");
        }
        else {
            out.print("GT");
        }
    }
}
