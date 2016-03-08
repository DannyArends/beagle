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
package blbutil;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.zip.GZIPOutputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;

/**
 * Class {@code FileUtil} contains static methods for working with files.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FileUtil {

    private FileUtil() {
        // private constructor prevents instantiation
    }

    /**
     * Returns a buffered {@code java.io.DataInputStream} reading from the
     * specified file.  If the input stream cannot be opened, an error message
     * will be printed and the java interpreter will exit.
     * @param file a file
     * @return a buffered {@code java.io.DataInputStream} reading from the
     * specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static DataInputStream dataInputStream(File file) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(
                    new FileInputStream(file)));
        } catch (FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return dis;
    }

    /**
     * Returns a buffered {@code java.io.DataOutputStream} writing to
     * the specified file.  Any existing file corresponding to the
     * {@code File} object will be deleted.   If the file cannot be opened,
     * an error message will be printed and the java interpreter will exit.
     * @param file a file
     * @return a buffered {@code java.io.DataOutputStream} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static DataOutputStream dataOutputStream(File file) {
        OutputStream dos = null;
        try {
            dos = new FileOutputStream(file);
        } catch (FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        DataOutputStream out = new DataOutputStream(
                new BufferedOutputStream(dos));
        return out;
    }

    /**
     * Returns a {@code java.io.PrintWriter} that writes
     * to standard out.
     *
     * @return a {@code java.io.PrintWriter} that writes
     * to standard out
     */
    public static PrintWriter stdOutPrintWriter() {
        return new PrintWriter(
                new BufferedOutputStream(System.out));
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to
     * the specified file.  The resulting file will be compressed using
     * the GZIP compression algorithm.  Any existing file corresponding
     * to the specified file will be deleted.  If the file
     * cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     * @param file a file
     * @return a {@code java.io.PrintWriter} writing to the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter gzipPrintWriter(File file) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(
                    new GZIPOutputStream(new FileOutputStream(file)));
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to
     * the specified file.  The resulting file will be compressed using
     * the BGZIP compression algorithm.  Any existing file corresponding
     * to the specified file will be deleted.  If the file
     * cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     *
     * @param file a file
     * @return a buffered {@code java.io.PrintWriter} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter bgzipPrintWriter(File file) {
        PrintWriter out = null;
        try {
            OutputStream fout = new FileOutputStream(file);
            out = new PrintWriter(new BlockCompressedOutputStream(
                    new BufferedOutputStream(fout), file));
        } catch (FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to
     * the specified file.  Any existing file corresponding
     * to the specified filename will be deleted.  If the file
     * cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     * @param file a file
     * @return a buffered {@code java.io.PrintWriter} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter printWriter(File file) {
        return printWriter(file, false);
    }

    /**
     * Returns a buffered {@code java.io.PrintWriter} writing to
     * the specified file. If {@code append == false}
     * any existing file corresponding to the specified file will be deleted.
     * If the file cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     *
     * @param file a file
     * @param append {@code true} if the data will be appended
     * to the end of any existing file
     * @return a buffered {@code java.io.PrintWriter} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter printWriter(File file, boolean append) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(
                    new BufferedWriter(new FileWriter(file, append)));
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns a non-buffered {@code java.io.PrintWriter} writing to
     * the specified file.
     * If {@code append == false} any existing file corresponding
     * to the specified file will be deleted. If the file cannot be opened,
     * an error message will be printed and the java interpreter will exit.
     *
     * @param file a file
     * @param append {@code true} if the data will be appended
     * to the end of any existing file
     * @return a non-buffered {@code java.io.PrintWriter} writing to
     * the specified file
     * @throws NullPointerException if {@code file == null}
     */
    public static PrintWriter nonBufferedPrintWriter(File file, boolean append) {
        boolean autoflush = true;
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileWriter(file, append), autoflush);
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return pw;
    }

    /**
     * Returns a temporary {@code File} that will be deleted when
     * the Java virtual machine exits.
     *
     * @param prefix the filename prefix.
     *
     * @return a {@code File} a new empty file.
     *
     * @throws IllegalArgumentException if {@code prefix} contains fewer than
     * three characters
     */
    public static File tempFile(String prefix) {
        File tempFile = null;
        try {
            tempFile = File.createTempFile(prefix, null);
            tempFile.deleteOnExit();
        } catch (IOException e) {
            Utilities.exit("Exception thrown by createTempFile: ", e);
        }
        return tempFile;
    }
}
