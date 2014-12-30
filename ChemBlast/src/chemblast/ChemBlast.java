/*
 // This Java source file is copyright (C) 2014 by Syed Asad Rahman. All rights
 // reserved. For further information, contact the author, Syed Asad Rahman, at
 // asad@ebi.ac.uk.
 //
 //
 // A copy of the GNU General Public License is provided in the file gpl.txt. You
 // may also obtain a copy of the GNU General Public License on the World Wide
 // Web at http://www.gnu.org/licenses/gpl.html.
 */
package chemblast;

import fingerprints.blast.Alignment;
import fingerprints.blast.BLAST;
import fingerprints.blast.BLASTFP;
import fingerprints.blast.DefaultAlignmentStats;
import fingerprints.blast.FingerprintDatabase;
import fingerprints.blast.FingerprintDatabaseGenerator;
import fingerprints.blast.FingerprintSequence;
import fingerprints.blast.Mol2Sequence;
import fingerprints.blast.ResultContainer;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

/**
 *
 * @author Syed Asad Rahman, e-mail: s9asad@gmail.com
 */
public class ChemBlast {

    final static SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    final static String newline = System.getProperty("line.separator");

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        ChemBlast chemBlast = new ChemBlast();
        List<String> list = Arrays.asList(args);
        if (list.contains("-formatDB") && args.length == 3) {
            if (list.contains("-DB")) {
                int index = list.indexOf("-DB") + 1;
                String dbName = list.get(index);
                File f = new File(dbName);
                if (!f.exists()) {
                    System.err.println("ERROR: The DB name " + f.getAbsolutePath() + " not found!");
                } else {
                    System.out.println("INFO: The DB name is " + f.getAbsolutePath());
                    chemBlast.formatDB(f);
                }
            } else {
                System.err.println("WARNING: -DB option missing! ");
                printHelp();
                System.exit(1);
            }
        } else if (list.contains("-searchDB") && args.length == 4) {
            int index = list.indexOf("-searchDB") + 1;
            String querySMI = null;
            if (list.contains("-querySMI")) {
                int smiIndex = list.indexOf("-querySMI") + 1;
                querySMI = list.get(smiIndex);
            } else {
                printHelp();
                System.err.println("WARNING: -querySMI option missing! ");
                System.exit(1);
            }
            if (querySMI == null) {
                System.err.println("ERROR: Missing input SMILES!");
                printHelp();
                System.exit(1);
            }
            String dbName = list.get(index);
            File databasefile = new File(dbName);
            if (!databasefile.exists()) {
                System.err.println("ERROR: The DB name " + databasefile.getAbsolutePath() + " not found!");
            } else {
                System.out.println("INFO: The DB name is " + databasefile.getAbsolutePath());
                String fileName = databasefile.getName().split("\\.")[0];
                File indexfile = new File(databasefile.getParentFile(), fileName + ".index");
                File formatDB = new File(databasefile.getParentFile(), fileName + ".format");
                if (!formatDB.exists() || !indexfile.exists()) {
                    chemBlast.formatDB(databasefile);
                } else {
                    System.out.println("INFO: Index file found " + indexfile.getAbsolutePath());
                    System.out.println("INFO: Formatted DB found " + formatDB.getAbsolutePath());
                }
                /*
                 Perform BLAST search
                 */
                try {
                    chemBlast.search(querySMI, indexfile, formatDB);
                } catch (InvalidSmilesException | IOException ex) {
                    Logger.getLogger(ChemBlast.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        } else {
            printHelp();
            System.exit(1);
        }
    }

    private void search(String querySmile, File indexfile, File formatDB)
            throws InvalidSmilesException, IOException {
        System.out.println("INFO: Searching DB...");
        int top_hits = 10;
        IAtomContainer parseSmiles = smilesParser.parseSmiles(querySmile);
        FingerprintSequence query = Mol2Sequence.convertMol2Sequence("Query", parseSmiles);
        long t1 = System.currentTimeMillis();

        FingerprintDatabase fpDB;
        fpDB = new FingerprintDatabase(formatDB, indexfile);
        System.out.println("INFO: Number of FP index: " + fpDB.getFingerprintCount());
        BLAST aligner = new BLASTFP(fpDB.getFingerprintCount());

        List<ResultContainer> container = new ArrayList<>();
        float percentage = 1.0f;
        for (long i = 0; i < fpDB.getFingerprintCount() * percentage; i++) {
            Alignment[] alignments = aligner.align(query, fpDB.getFingerprintSequence(i));
            ResultContainer alignmentContainer;
            DefaultAlignmentStats defaultAlignmentStats;
            // Print the alignments
            for (Alignment alignment : alignments) {
                defaultAlignmentStats = new DefaultAlignmentStats(alignment.getSubjectLength());
                alignmentContainer
                        = new ResultContainer(defaultAlignmentStats, alignment, query, fpDB.getFingerprintSequence(i));
                container.add(alignmentContainer);
            }
        }
        long t2 = System.currentTimeMillis();
        System.out.printf("Running Time: %d msec%n", t2 - t1);
        fpDB.close();

        /*
         Sort the results
         */
        Collections.sort(container);

        if (container.size() > top_hits) {
            System.out.println("INFO: Top " + top_hits + " Hits");
        } else if (container.size() <= top_hits && container.size() > 0) {
            System.out.println("INFO: Top Hits");
        } else {
            System.out.println("INFO: No significant hits found!");
        }
        /*
         Print sorted results
         */
        for (int i = 0; i < container.size(); i++) {
            ResultContainer alignmentContainer = container.get(i);
            System.out.println("Query " + alignmentContainer.getQuery().description()
                    + ", Subject " + alignmentContainer.getSubject().description()
                    + ", bit-Score "
                    + alignmentContainer.getBitScore()
                    + ", rawScore " + alignmentContainer.getRawScore()
                    + ", query-Coverage "
                    + alignmentContainer.getPercentageQueryCoverage()
                    + "%, eValue "
                    + alignmentContainer.getEValue()
                    + ", Identity "
                    + alignmentContainer.getPercentageIdentities()
                    + "%, Positives "
                    + alignmentContainer.getPercentagePositives()
                    + "%");
            //alignmentContainer.printDetails(System.out);
            if (i == top_hits) {
                break;
            }
        }
    }

    private static void printHelp() {
        System.out.println("HELP!");
        System.out.println("1) Formatting the input database. Each line in the input file contains identifier for the molecule/drug and its smiles (Name\tSMILES).");
        System.out.println("\tjava -jar -formatDB -DB smiles.txt");
        System.out.println("2) Searching the input database. Input smiles file and the query smiles. This will return max. of 10 hits");
        System.out.println("\tjava -jar -searchDB smiles.txt -querySMI \"N1=C(C#N)N=C2C(=C1NCC(CC)C)NCN2CC3=CC=CC=C3\"");
    }

    public ChemBlast() {
    }
    /*
     Format the input DB- a) generate SMILES and b) transform it into tuples
     */

    void formatDB(File databasefile) {
        System.out.println("INFO: Formatting input database....");
        try {
            String fileName = databasefile.getName().split("\\.")[0];
            File indexfile = new File(databasefile.getParentFile(), fileName + ".index");
            File formatDB = new File(databasefile.getParentFile(), fileName + ".format");
            String numberOfSeq = "-1";
            FingerprintDatabaseGenerator fingerprintDatabaseGenerator
                    = new FingerprintDatabaseGenerator(databasefile, indexfile, formatDB, numberOfSeq);
            System.out.println("INFO: Index file created " + indexfile.getAbsolutePath());
            System.out.println("INFO: DB formatted " + formatDB.getAbsolutePath());
            System.out.println("INFO: Done");
        } catch (IOException ex) {
            Logger.getLogger(ChemBlast.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
