package ru.bmstu.bioinformatics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class BLOSUM62 implements ScoringFunction {
    private final String BLOSUM_62_MATRIX_PATH = "BLOSUM62.dat", EMPTY = "", SPACE_REGEX = "\\s+";
    private final int NULL_DELTA = 0, SINGLE_DELTA = 1, START_CHARACTER_ORDER = 0,
            FIRST_ALPHABET_ELEMENT = 0, FIRST_INDEX = 1, CHARACTER_INDEX = 0, EXIT_CODE_WITH_ERROR = 1;

    private int indel;

    private HashMap<Character, Integer> orderInMatrix = new HashMap<>();
    private ArrayList<ArrayList<Integer>> matrix = new ArrayList<>();

    private void initializeMatrix() {
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(BLOSUM_62_MATRIX_PATH));

            String line = bufferedReader.readLine();
            String[] alphabet = line.split(SPACE_REGEX);
            int delta = NULL_DELTA;
            if (alphabet[FIRST_ALPHABET_ELEMENT].equals(EMPTY))
                delta = SINGLE_DELTA;

            for (int order = START_CHARACTER_ORDER + delta; order < alphabet.length; order++) {
                orderInMatrix.put(alphabet[order].charAt(CHARACTER_INDEX), order - delta);
            }
            while ((line = bufferedReader.readLine()) != null) {
                String[] splitLine = line.split(SPACE_REGEX);
                ArrayList<Integer> scores = new ArrayList<>();
                for (int index = FIRST_INDEX; index < splitLine.length; index++) {
                    scores.add(Integer.parseInt(splitLine[index]));
                }
                matrix.add(scores);
            }

            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(EXIT_CODE_WITH_ERROR);
        }
    }

    public BLOSUM62(int indel) {
        this.indel = indel;
        initializeMatrix();
    }

    public int getIndel() {
        return indel;
    }

    public int getScore(char firstChar, char secondChar) {
        int firstOrder = orderInMatrix.get(firstChar),
                secondOrder = orderInMatrix.get(secondChar);
        return matrix.get(firstOrder).get(secondOrder);
    }
}
