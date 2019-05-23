#!/usr/bin/env python
# coding: utf-8

# In[12]:


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
    plt.axis([0, (len(cover) // 1000 - 1) * 1000, 0, max(avcover)+10])
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(filename)
    plt.show()


# In[13]:


filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW3/Data/pb_aln.bam'

print(showCover(filename))


# In[14]:


filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW3/Data/ont_aln.bam'

print(showCover(filename))


# In[15]:


import pysam
import matplotlib.pyplot as plt
# import Bio
from Bio import SeqIO
import numpy as np
import math

import math
import pandas as pd


def count_subst_freqs(filename, dataset):
    freqs_matrix = {}
    for ref in ('ACTG-'):
        for alt in ('ACTG-'):
            if alt != ref:
                freqs_matrix[(alt, ref)] = 0
                
    bamfile = pysam.AlignmentFile(filename, 'rb')
    for read in bamfile.fetch():
        try:
            aln_pairs = read.get_aligned_pairs(with_seq=True)
        except:
            continue
        
        # remove skippings
        s1 = -1
        s2 = len(aln_pairs)
        while aln_pairs[s1 + 1][1] is None:
            s1 += 1
        while aln_pairs[s2 - 1][1] is None:
            s2 -= 1
        aln_pairs = aln_pairs[s1 + 1 : s2]
            
            
        read_seq = read.seq
        for p in aln_pairs:
            if p[0] is None:
                freqs_matrix[('-', p[2])] += 1
            elif p[1] is None and read_seq[p[0]] != 'N':
                freqs_matrix[(read_seq[p[0]], '-')] += 1
            elif p[1] is not None and p[2].lower() == p[2] and read_seq[p[0]] != 'N':
                freqs_matrix[(read_seq[p[0]], p[2].upper())] += 1
                
    print('File: {}\n'.format(dataset))
    print('   A\t\tC\t\tT\t\tG\t\t-')
    for alt in 'ACTG-':
        print(alt, end='  ')
        for ref in 'ACTG-':
            if alt == ref:
                print('--\t\t', end='')
            else:
                print(freqs_matrix[(alt, ref)], 
                      end='\t' + ('\t' if len(str(freqs_matrix[(alt, ref)])) < 5 or ref != 'A' else ''))
        print()


# In[22]:


# filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2/Data/MG1655-K12.first10K.fasta'
sam_filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW3/Data/ont_aln_with_MD.bam'


count_subst_freqs(sam_filename, 'ont_aln.bam')


# In[21]:


sam_filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW3/Data/pb_aln_with_MD.bam'


count_subst_freqs(sam_filename, 'pb_aln')


# In[23]:


def indel_stats(filename, dataset):
    length = []
    ins = []
    dels = []
    
    with open(filename) as f:
        for line in f:
            line = line.replace('\n', '').split('\t')
            length.append(int(line[0]))
            ins.append(int(line[1]))
            dels.append(int(line[2]))
    
    mean_ins = sum([i * l for (i, l) in zip(ins, length)]) / sum(ins)
    mean_dels = sum([d * l for (d, l) in zip(dels, length)]) / sum(dels)
    
    print('Data: {}\n'
          'Mean insertion len: {:.3f}\n'
          'Mean deletion len: {:.3f}\n'
          .format(dataset, mean_ins, mean_dels))
    
    plt.figure(figsize=(7, 12))
    plt.subplot(211)
    plt.plot(length, ins)
    plt.xlabel('length')
    plt.ylabel('number')
    plt.title('Insertion')
    
    plt.subplot(212)
    plt.plot(length, dels)
    plt.xlabel('length')
    plt.ylabel('number')
    plt.title('Deletion')


# In[24]:


indel_stats('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW3/Data/pb_indel_distr', 'pb')


# In[25]:


indel_stats('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW3/Data/ont_indel_distr', 'ont')


# In[ ]:




