import pysam
import matplotlib.pyplot as plt
%matplotlib inline
from Bio import SeqIO
import pandas as pd


def calculateStatistic(record1, record2, ref):
    if record1.is_unmapped or record2.is_unmapped:
        return (0, 0, 0)
    if record1.reference_start != record2.reference_start:
        return (0, 0, 0)

    FP = 0
    FN = 0
    TN = 0

    cigar1 = record1.cigartuples
    cigar2 = record2.cigartuples
    rpos = record1.reference_start
    seq1 = record1.query_sequence
    seq2 = record2.query_sequence
    spos1 = 0
    spos2 = 0
    ci1 = 0
    ci2 = 0
    cnt1 = 0
    cnt2 = 0
    while (ci1 < len(cigar1)) or (ci2 < len(cigar2)):
        if (ci1 < len(cigar1) and (cigar1[ci1][0] == 0 or cigar1[ci1][0] == 7 or cigar1[ci1][0] == 8)) and (
                ci2 < len(cigar2) and (cigar2[ci2][0] == 0 or cigar2[ci2][0] == 7 or cigar2[ci2][0] == 8)):
            if (seq1[spos1] == ref[rpos]) and (seq2[spos2] != ref[rpos]):
                FN += 1
            elif (seq1[spos1] != ref[rpos]) and (seq2[spos2] == ref[rpos]):
                TN += 1
            elif (seq1[spos1] != ref[rpos]) and (seq2[spos2] != ref[rpos]):
                FP += 1
            cnt1 += 1
            cnt2 += 1
            rpos += 1
            spos1 += 1
            spos2 += 1
        elif (ci1 < len(cigar1) and (cigar1[ci1][0] == 1 or cigar1[ci1][0] == 4)) and (
                (ci2 < len(cigar2)) and (cigar2[ci2][0] == 1 or cigar2[ci2][0] == 4)):
            FP += 1
            cnt1 += 1
            cnt2 += 1
            spos1 += 1
            spos2 += 1
        elif (ci1 < len(cigar1) and (cigar1[ci1][0] == 1 or cigar1[ci1][0] == 4)):
            TN += 1
            cnt1 += 1
            spos1 += 1
        elif (ci2 < len(cigar2) and (cigar2[ci2][0] == 1 or cigar2[ci2][0] == 4)):
            FN += 1
            cnt2 += 1
            spos2 += 1
        else:
            if (ci1 < len(cigar1) and (cigar1[ci1][0] == 0 or cigar1[ci1][0] == 7 or cigar1[ci1][0] == 8)):
                spos1 += 1
            if (ci2 < len(cigar2) and (cigar2[ci2][0] == 0 or cigar2[ci2][0] == 7 or cigar2[ci2][0] == 8)):
                spos2 += 1
            if (ci1 < len(cigar1)):
                cnt1 += 1
            if (ci2 < len(cigar2)):
                cnt2 += 1
            rpos += 1

        if (ci1 < len(cigar1) and cnt1 == cigar1[ci1][1]):
            cnt1 = 0
            ci1 += 1
        if (ci2 < len(cigar2) and cnt2 == cigar2[ci2][1]):
            cnt2 = 0
            ci2 += 1
    return (FP, FN, TN)


def main(refFileName, oldFileName, newFileName):
    ref = ""
    for rec in SeqIO.parse(refFileName, "fasta"):
        ref = rec.seq
    oldSamFile = pysam.AlignmentFile(oldFileName, "r")
    newSamFile = pysam.AlignmentFile(newFileName, "r")

    FP = 0
    FN = 0
    TN = 0
    iter1 = oldSamFile.fetch()
    cnt = 0
    for record2 in newSamFile.fetch():
        if (cnt % 10000 == 0):
            print(cnt, record2.query_name)
        record1 = next(iter1)
        while (record1.query_name != record2.query_name):
            record1 = next(iter1)
        cFP, cFN, cTP = calculateStatistic(record1, record2, ref)
        FP += cFP
        FN += cFN
        TN += cTP
        cnt += 1