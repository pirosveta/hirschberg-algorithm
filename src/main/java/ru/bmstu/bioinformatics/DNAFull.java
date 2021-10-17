package ru.bmstu.bioinformatics;

public class DNAFull implements ScoringFunction {
    private int indel;
    private int match;
    private int mismatch;

    public DNAFull(int indel, int match, int mismatch) {
        this.indel = indel;
        this.match = match;
        this.mismatch = mismatch;
    }

    public int getIndel() {
        return indel;
    }

    public int getMatch() {
        return match;
    }

    public int getMismatch() {
        return mismatch;
    }
}
