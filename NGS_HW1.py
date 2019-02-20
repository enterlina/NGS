from Bio import SeqIO
import matplotlib.pyplot as plt



def QualityRead(cntAT, cntGC, cntN, bucket_cnt):
    if (cntAT + cntGC > bucket_cnt):
        return True
    return False


def GCDistribution(file_name, bucket_cnt):
    ans = [0] * (bucket_cnt)
    cnt = 0
    for rec in SeqIO.parse(file_name, "fastq"):
        cnt += 1
        cntAT = 0
        cntGC = 0
        cntN = 0
        reclen = len(rec.seq)
        for i in range(reclen):
            quality_trash_hold = 30
            if (rec.letter_annotations["phred_quality"][i] > quality_trash_hold):
                if (rec.seq[i] == 'A' or rec.seq[i] == 'T'):
                    cntAT += 1
                elif (rec.seq[i] == 'G' or rec.seq[i] == 'C'):
                    cntGC += 1
                elif (rec.seq[i] == 'N'):
                    cntN += 1
            else:
                cntN += 1

        if QualityRead(cntAT, cntGC, cntN, bucket_cnt):
            if cntAT != 0:
                ans[cntGC * bucket_cnt // (cntGC + cntAT)] += 1
    return ans


# import matplotlib.pyplot as plt
# # % matplotlib
# # inline


def resGCDistribution(file_name, bcnt):
    delta = 100 // bcnt
    xs = [(i + (i + delta)) / 2 for i in range(0, 100, delta)]
    ys = GCDistribution(file_name, bcnt)

    print(xs)
    print(ys)

    plt.plot(xs, ys, 'r+')
    plt.axis([0, 100, 0, max(ys)+2])
    plt.xlabel('%GC')
    plt.ylabel('Number of reads')
    plt.title(file_name)
    plt.show()

resGCDistribution("/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW1/test.fastq", 20)

# resGCDistribution("/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW1/ecoli_mda_lane1_right.fastq", 20)
