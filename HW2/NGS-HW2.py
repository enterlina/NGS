import pysam
import matplotlib.pyplot as plt
# import Bio
from Bio import SeqIO
import numpy as np
import math

import math
import pandas as pd


def Cover(filename):
    samfile = pysam.AlignmentFile(filename, "r")
    reflen = samfile.header.get("SQ")[0]['LN']

    cover_len = [0] * (reflen + 1)

    for record in samfile.fetch():
        if not record.is_unmapped:
            rpos = record.reference_start
            rlen = record.reference_length
            cover_len[rpos] += 1
            cover_len[rpos + rlen] -= 1

    for i in range(1, reflen):
        cover_len[i] += cover_len[i - 1]
    return cover_len


filename = '/Users/alena_paliakova/Desktop/HW2/MG1655-K12_aligned.sam'


# filename = '/Users/alena_paliakova/Desktop/HW2/MG1655-K12_sorted.bam'
# print(Cover(filename))


def showCover(filename):
    cover = Cover(filename)
    print("Average coverage: " + str(sum(cover) / len(cover)))
    print("Coverage rate: " + str((len(cover) - cover.count(0)) / len(cover)))
    avcover = [0] * (len(cover) // 1000 + 1)
    x = [i for i in range(0, len(cover), 1000)]
    for i in range(len(cover)):
        avcover[i // 1000] += cover[i]
    for i in range(len(avcover) - 1):
        avcover[i] /= 1000
    avcover = avcover[:-1]
    x = x[:len(avcover)]

    plt.plot(x, avcover)
    plt.axis([0, (len(cover) // 1000 - 1) * 1000, 0, max(avcover) + 200])
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(filename)
    plt.show()


# print(showCover(filename))


def getInsertSize(filename):
    samfile = pysam.AlignmentFile(filename, "r")
    insertSize = []

    for record in samfile.fetch():
        if ((not record.is_unmapped) and record.is_paired and record.is_read1):
            if (record.template_length > 0):
                insertSize.append(record.template_length)
            # rpos = record.reference_start
            # rlen = record.reference_length

            # npos = record.next_reference_start
            # nlen = record.next_reference_length
            # if (npos > (rpos + rlen)):
            #    insertSize.append(npos + nlen - rpos)
    return insertSize


def showInsertSize(filename):
    insertSize = getInsertSize(filename)
    npis = np.array(insertSize)
    print("Average insert size: " + str(np.mean(npis)))
    print("SD insert size: " + str(np.std(npis)))
    insertSize.sort()
    bbg = 0
    cntElem = int(0.95 * len(insertSize))
    for bg in range(int(len(insertSize) - 0.95 * len(insertSize))):
        if ((insertSize[bg + cntElem] - insertSize[bg]) < insertSize[bbg + cntElem] - insertSize[bbg]):
            bbg = bg
    print("95% seq: " + str(insertSize[bbg]) + "-" + str(insertSize[bbg + cntElem]))

    for ig in range(0, 1):
        npis = np.array(insertSize)
        lst = int(len(insertSize))
        mn = insertSize[len(insertSize) // 2] - np.std(npis) ** 0.5
        mx = insertSize[len(insertSize) // 2] + np.std(npis) ** 0.5
        print(mn, mx)
        bg = 0
        for i in range(int(len(insertSize))):
            if (insertSize[i] < mx):
                lst = i
            if (insertSize[i] < mn):
                bg = i

        insertSize = insertSize[bg:lst]

    num_bins = 50
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(insertSize, num_bins)

    ax.set_xlabel('Insert size')
    ax.set_ylabel('Reads number')
    ax.set_title(filename)
    # ax.set_xlim([insertSize[bbg] - 20, insertSize[bbg + cntElem] + 20])

    fig.tight_layout()
    # mpl.rcParams['figure.figsize'] = (100,100)
    plt.show()


# print ('Insert sine', getInsertSize(filename))
# print (showInsertSize(filename))


def getLetId(c, isrv):
    # print(c)
    let = ["ACGTN", "TGCAN"]
    for i in range(5):
        if (let[isrv][i] == c):
            return i


def getReplacementMatrix(reffilename, samfilename):
    for rec in SeqIO.parse(reffilename, "fasta"):
        ref = rec.seq
    # ref = rec.seq
    samfile = pysam.AlignmentFile(samfilename, "r")
    res = [[0] * 5 for i in range(5)]
    let = "ACGTN"

    for record in samfile.fetch():
        if record.is_unmapped:
            continue
        cigar = record.cigartuples
        rpos = record.reference_start
        seq = record.query_sequence
        spos = 0
        for i in range(len(cigar)):
            if (cigar[i][0] == 0 or cigar[i][0] == 7 or cigar[i][0] == 8):
                cntl = cigar[i][1]
                isrv = 0  # record.is_reverse
                for j in range(cntl):
                    res[getLetId(seq[spos + j], isrv)][getLetId(ref[rpos + j], 0)] += 1
                rpos += cigar[i][1]
                spos += cigar[i][1]
            elif (cigar[i][0] == 1 or cigar[i][0] == 4):
                spos += cigar[i][1]
            elif (cigar[i][0] == 2 or cigar[i][0] == 3):
                rpos += cigar[i][1]
    samfile.close()
    res = res[:4]
    for i in range(4):
        res[i] = res[i][:4]
    return res


filename = "/Users/alena_paliakova/Desktop/HW2/MG1655-K12.first10K.fasta"
samfile = "MG1655-K12_aligned.sam"
res = getReplacementMatrix(filename, samfile)

df = pd.DataFrame(data=res, columns=['A', 'C', 'G', 'T'], index=['A', 'C', 'G', 'T'])
print()
