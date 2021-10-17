package ru.bmstu.bioinformatics;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class Main {
    @Parameter(names = "-i", description = "Paths to input sequences", arity = 2, required = true)
    private static List<String> inputPaths = new ArrayList<>();

    @Parameter(names = "-a", description = "Alphabet")
    private static String alphabet = "Default";

    @Parameter(names = {"-g", "--gap"}, description = "Penalty for the gap")
    private static String gapPenalty = "";

    @Parameter(names = "-o", description = "Path to output file")
    private static String outputPath = "";

    private static final String DNA_FULL = "DNAFull", BLOSUM_62 = "BLOSUM62", DEFAULT = "Default", EMPTY = "";
    private static final int MATCH_DNA_FULL = 5, MISMATCH_DNA_FULL = -4,
            INDEL_DEFAULT = -2, MATCH_DEFAULT = 1, MISMATCH_DEFAULT = -1,
            NUMBER_OF_SEQUENCES = 2, FIRST_SEQUENCE_INDEX = 0, SECOND_SEQUENCE_INDEX = 1,
            EXIT_CODE_WITH_ERROR = 1;

    public static ArrayList<String> readFile() {
        ArrayList<String> sequences = new ArrayList<>();

        for (String sequencePath : inputPaths) {
            try {
                BufferedReader bufferedReader = new BufferedReader(new FileReader(sequencePath));
                StringBuilder stringBuilder = new StringBuilder();
                String line;

                bufferedReader.readLine();
                while ((line = bufferedReader.readLine()) != null) {
                    stringBuilder.append(line);
                }
                sequences.add(stringBuilder.toString());

                bufferedReader.close();
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(EXIT_CODE_WITH_ERROR);
            }
        }

        return sequences;
    }

    private static ScoringFunction defineScoringFunction() {
        ScoringFunction scoringFunction;

        if (alphabet.equals(DNA_FULL)) {
            scoringFunction = new DNAFull(Integer.parseInt(gapPenalty), MATCH_DNA_FULL, MISMATCH_DNA_FULL);
        } else if (alphabet.equals(BLOSUM_62)) {
            scoringFunction = new BLOSUM62(Integer.parseInt(gapPenalty));
        } else scoringFunction = new Default(INDEL_DEFAULT, MATCH_DEFAULT, MISMATCH_DEFAULT);

        return scoringFunction;
    }

    public static void main(String[] args) {
        try {
            Main main = new Main();
            JCommander jCommander = new JCommander(main);
            jCommander.parse(args);

            if (inputPaths.size() != NUMBER_OF_SEQUENCES
                    || !(alphabet.equals(BLOSUM_62) || alphabet.equals(DNA_FULL) || alphabet.equals(DEFAULT))
                    || (!alphabet.equals(DEFAULT) && gapPenalty.equals(EMPTY))) {
                jCommander.usage();
                return;
            }

            ArrayList<String> sequences = readFile();
            PairAlignment pairAlignment = new PairAlignment(
                    sequences.get(FIRST_SEQUENCE_INDEX),
                    sequences.get(SECOND_SEQUENCE_INDEX),
                    defineScoringFunction());

            PrintStream console = System.out;
            if (!outputPath.equals(EMPTY)) {
                File file = new File(outputPath);
                FileOutputStream fos = new FileOutputStream(file);
                PrintStream ps = new PrintStream(fos);
                System.setOut(ps);
            }
            System.out.println(pairAlignment.toString());
            System.setOut(console);

        } catch (ParameterException | FileNotFoundException e) {
            e.printStackTrace();
        }

    }
}