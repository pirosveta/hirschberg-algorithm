package ru.bmstu.bioinformatics;

import java.math.BigInteger;
import java.util.*;

public class PairAlignment {
    private final String OUTPUT_SEQUENCE_1 = "Seq1: ", OUTPUT_SEQUENCE_2 = "Seq2: ", OUTPUT_SCORE = "Score: ",
            NEXT_LINE = "\n";
    private final char GAP = '_';
    private final int ZERO = 0, GAP_INDEX = 0, ONE = 1, FIRST_INDEX = 1, NUMBER_OF_SYMBOLS_IN_LINE = 50;

    private TreeSet<int[]> path = new TreeSet<>(new Comparator<int[]>() {
        @Override
        public int compare(int[] o1, int[] o2) {
            if (o1[0] < o2[0] && o1[1] < o2[1]) {
                return -1;
            } else if (o2[0] < o1[0] && o2[1] < o1[1]) {
                return 1;
            } else if (o1[0] == o2[0] && o1[1] == o2[1]) {
                return 0;
            } else if (o1[0] == o2[0]) {
                return (o1[1] < o2[1] ? -1 : 1);
            } else return (o1[0] < o2[0] ? -1 : 1);
        }
    });
    private ScoringFunction scoringFunction;

    private StringBuilder firstAlignSequence = new StringBuilder(),
            secondAlignSequence = new StringBuilder();
    private Integer score = null;

    private int getIndel() {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getIndel();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getIndel();
        } else return ((BLOSUM62) scoringFunction).getIndel();
    }

    private int getMatch(char firstChar, char secondChar) {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getMatch();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getMatch();
        } else return ((BLOSUM62) scoringFunction).getScore(firstChar, secondChar);
    }

    private int getMismatch(char firstChar, char secondChar) {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getMismatch();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getMismatch();
        } else return ((BLOSUM62) scoringFunction).getScore(firstChar, secondChar);
    }

    private boolean isDegenerateCase(StringBuilder firstSequence, StringBuilder secondSequence) {
        return firstSequence.length() < 2 || secondSequence.length() < 2;
    }

    private int getMaximumScoreIndex(ArrayList<Integer> firstPart, ArrayList<Integer> secondPart) {
        int maximumScore = Integer.MIN_VALUE, indexScore = 0, partSize = firstPart.size() - 1;

        int counter = secondPart.size() - 1;
        for (int firstElement : firstPart) {
            int secondElement = secondPart.get(counter);
            if (firstElement + secondElement > maximumScore) {
                maximumScore = firstElement + secondElement;
                indexScore = partSize - counter;
            }
            counter--;
        }

        if (score == null) {
            score = maximumScore;
        }
        return indexScore - 1;
    }

    private ArrayList<Integer> findScore(String firstSequence, String secondSequence) {
        int bufferSize = secondSequence.length() + 1, cycleSize = firstSequence.length() + 1,
                columnIndex = 1, maximumScore, nextScore,
                diagElement, leftElement, upElement;
        char firstChar, secondChar;
        ArrayList<Integer> firstBuffer = new ArrayList<>(bufferSize),
                secondBuffer = new ArrayList<>(bufferSize);

        for (int index = 0; index < bufferSize; index++) {
            firstBuffer.add(getIndel() * index);
        }
        while (columnIndex < cycleSize) {
            diagElement = firstBuffer.get(0);
            upElement = getIndel() * columnIndex;
            secondBuffer.add(upElement);

            for (int lineIndex = 1; lineIndex < bufferSize; lineIndex++) {
                leftElement = firstBuffer.get(lineIndex);
                firstChar = firstSequence.charAt(columnIndex - 1);
                secondChar = secondSequence.charAt(lineIndex - 1);
                boolean isMatch = firstChar == secondChar;

                maximumScore = diagElement
                        + (isMatch ? getMatch(firstChar, secondChar) : getMismatch(firstChar, secondChar));
                nextScore = leftElement + getIndel();
                if (nextScore > maximumScore) maximumScore = nextScore;
                nextScore = upElement + getIndel();
                if (nextScore > maximumScore) maximumScore = nextScore;

                secondBuffer.add(maximumScore);

                upElement = maximumScore;
                diagElement = leftElement;
            }

            firstBuffer.clear();
            firstBuffer.addAll(secondBuffer);
            secondBuffer.clear();

            columnIndex++;
        }

        return firstBuffer;
    }

    private int[] getIndexes(StringBuilder splitString, StringBuilder wholeString) {
        int splitStringLength = splitString.length(),
                splitStringIndex = (splitStringLength % 2 == 0)
                        ? (splitStringLength / 2 - 1)
                        : (splitStringLength / 2),
                wholeStringIndex;

        ArrayList<Integer> firstPartScore = findScore(
                splitString.substring(0, splitStringIndex + 1),
                wholeString.toString());
        ArrayList<Integer> secondPartScore = findScore(
                splitString.reverse().substring(0, splitStringLength - (splitStringIndex + 1)),
                wholeString.reverse().toString());
        wholeStringIndex = getMaximumScoreIndex(firstPartScore, secondPartScore);

        return new int[]{splitStringIndex, wholeStringIndex};
    }

    private void fillPath(String firstSequenceString, String secondSequenceString, int lineDelta, int columnDelta) {
        StringBuilder firstSequence = new StringBuilder(firstSequenceString),
                secondSequence = new StringBuilder(secondSequenceString);
        int firstLength = firstSequence.length(), secondLength = secondSequence.length(),
                lineIndex, columnIndex;
        notDegenerateCase:
        if (!isDegenerateCase(firstSequence, secondSequence)) {
            if (firstLength >= secondLength) {
                int[] indexes = getIndexes(firstSequence, secondSequence);
                columnIndex = indexes[0];
                lineIndex = indexes[1];
            } else {
                int[] indexes = getIndexes(secondSequence, firstSequence);
                lineIndex = indexes[0];
                columnIndex = indexes[1];
            }

            path.add(new int[]{lineIndex + lineDelta, columnIndex + columnDelta});

            if (lineIndex < 1 && columnIndex < 1 && firstLength == 2 && secondLength == 2) {
                path.add(new int[]{lineIndex + 1 + lineDelta, columnIndex + 1 + columnDelta});
                break notDegenerateCase;
            }
            if (lineIndex == secondLength - 1 && columnIndex == firstLength - 1
                    && firstLength == 2 && secondLength == 2) {
                path.add(new int[]{lineIndex - 1 + lineDelta, columnIndex - 1 + columnDelta});
                break notDegenerateCase;
            }

            fillPath(firstSequenceString.substring(0, columnIndex + 1),
                    secondSequenceString.substring(0, lineIndex + 1),
                    lineDelta, columnDelta);
            if (lineIndex == -1) lineIndex = 0;
            if (columnIndex == -1) columnIndex = 0;
            fillPath(firstSequenceString.substring(columnIndex),
                    secondSequenceString.substring(lineIndex),
                    lineDelta + lineIndex, columnDelta + columnIndex);
        }
    }

    private void align(String firstSequence, String secondSequence) {
        int[] p = path.first();
        int lineIndex = p[0], prevLineIndex = lineIndex,
                columnIndex = p[1], prevColumnIndex = columnIndex;

        if (lineIndex != -1 && columnIndex != -1) {
            path.add(new int[]{lineIndex - 1, columnIndex - 1});
            lineIndex--;
            prevLineIndex = lineIndex;
            columnIndex--;
            prevColumnIndex = columnIndex;
        }
        if (lineIndex != -1) {
            char[] addition = new char[lineIndex + 1];
            Arrays.fill(addition, GAP);
            firstAlignSequence.append(addition);
            secondAlignSequence.append(new StringBuilder(secondSequence).substring(0, lineIndex + 1));
        }
        if (columnIndex != -1) {
            char[] addition = new char[columnIndex + 1];
            Arrays.fill(addition, GAP);
            firstAlignSequence.append(new StringBuilder(firstSequence).substring(0, columnIndex + 1));
            secondAlignSequence.append(addition);
        }

        Iterator<int[]> pathIterator = path.iterator();
        pathIterator.next();
        while (pathIterator.hasNext()) {
            p = pathIterator.next();
            lineIndex = p[0];
            columnIndex = p[1];

            if (lineIndex > prevLineIndex) {
                secondAlignSequence.append(secondSequence.charAt(lineIndex));
            } else {
                secondAlignSequence.append(GAP);
            }
            if (columnIndex > prevColumnIndex) {
                firstAlignSequence.append(firstSequence.charAt(columnIndex));
            } else {
                firstAlignSequence.append(GAP);
            }

            prevLineIndex = lineIndex;
            prevColumnIndex = columnIndex;
        }

        if (lineIndex < secondSequence.length() - 1) {
            char[] addition = new char[secondSequence.length() - 1 - lineIndex];
            Arrays.fill(addition, GAP);
            firstAlignSequence.append(addition);
            secondAlignSequence.append(secondSequence.substring(lineIndex + 1));
        }
        if (columnIndex < firstSequence.length() - 1) {
            char[] addition = new char[firstSequence.length() - 1 - columnIndex];
            Arrays.fill(addition, GAP);
            firstAlignSequence.append(firstSequence.substring(columnIndex + 1));
            secondAlignSequence.append(addition);
        }
    }

    public PairAlignment(String firstSequence, String secondSequence,
                         ScoringFunction scoringFunction) {
        this.scoringFunction = scoringFunction;

        fillPath(firstSequence, secondSequence, 0, 0);
        align(firstSequence, secondSequence);
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        int i = FIRST_INDEX;
        for (; i < Math.ceil((double) firstAlignSequence.length() / NUMBER_OF_SYMBOLS_IN_LINE); i++) {
            stringBuilder.append(OUTPUT_SEQUENCE_1);
            stringBuilder.append(firstAlignSequence.substring(
                    NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE),
                    NUMBER_OF_SYMBOLS_IN_LINE * i));
            stringBuilder.append(NEXT_LINE);

            stringBuilder.append(OUTPUT_SEQUENCE_2);
            stringBuilder.append(secondAlignSequence.substring(
                    NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE),
                    NUMBER_OF_SYMBOLS_IN_LINE * i));
            stringBuilder.append(NEXT_LINE);
            stringBuilder.append(NEXT_LINE);
        }

        stringBuilder.append(OUTPUT_SEQUENCE_1);
        stringBuilder.append(firstAlignSequence.substring(NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE)));
        stringBuilder.append(NEXT_LINE);

        stringBuilder.append(OUTPUT_SEQUENCE_2);
        stringBuilder.append(secondAlignSequence.substring(NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE)));
        stringBuilder.append(NEXT_LINE);
        stringBuilder.append(NEXT_LINE);

        stringBuilder.append(OUTPUT_SCORE);
        stringBuilder.append(score);

        return stringBuilder.toString();
    }
}
