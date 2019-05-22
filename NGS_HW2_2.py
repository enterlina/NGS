
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




# In[ ]:





# In[4]:


# filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2/Data/MG1655-K12.first10K.fasta'
sam_filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/b22_aln.bam'


count_subst_freqs(sam_filename, 'b22_aln.bam')



# In[5]:



sam_filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/c24_aln.bam'


count_subst_freqs(sam_filename, 'c24_aln.bam')


sam_filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/test_aln.bam'


count_subst_freqs(sam_filename, 'test_aln.bam')




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
    plt.axis([0, (len(cover) // 1000 - 1) * 1000, 0, max(avcover) + 200])
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(filename)
    plt.show()



# In[50]:


filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/b22_aln.bam'

print(showCover(filename))


# In[51]:


filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/c24_aln.bam'
print(showCover(filename))


# In[46]:


def count_error_percentage(filename, dataset):
    errors = []
    with open(filename) as f:
        for line in f:
            line = line.split(' ')
            # in other case read didn't aligned
            if line[1].replace('\n', '').isdigit():
                errors.append(int(line[1]) / int(line[0]))
    
    avg = sum(errors) / len(errors)
    
    print("Data: {}\n"
          "Mean errors %: {:.3f}%\n"
          .format(dataset, avg * 100))


# In[13]:



count_error_percentage('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/b22_errors', 'b22')


# In[19]:


count_error_percentage('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/c24_errors', 'c24')


# In[23]:


def statics_indels(filename, dataset):
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


# In[25]:


statics_indels('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/b22_indel_distr', 'b22')


# In[26]:


statics_indels('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/c24_indel_distr', 'c24')


# In[29]:


def homo_quality_indel(filename, dataset):
    hp_indel = {}
    for k in range(3, 11):
        hp_indel[k - 3] = {}
        for g in range(-25, 25):
            hp_indel[k - 3][g] = 0
        
    bamfile = pysam.AlignmentFile(filename, 'rb')
    curr_len = 1
    curr_seg_len = 0
    seg_len_pool = []
    curr_nucl = 'R'
    for read in bamfile.fetch():
        if curr_len > 2  and curr_len < 11:
            hp_indel[curr_len - 3][curr_seg_len - curr_len] += 1
            for l in seg_len_pool:
                hp_indel[curr_len - 3][l - curr_len] += 1
            
        curr_len = 1 
        curr_seg_len = 0
        curr_nucl = 'R'
        seg_len_pool = []
        try:
            aln_pairs = read.get_aligned_pairs(with_seq=True)
        except:
            continue
        
        for p in read.get_aligned_pairs(with_seq=True):
            if p[1] is None:
                if read.seq[p[0]] == curr_nucl:
                    curr_seg_len += 1
                else:
                    if curr_seg_len != 0:
                        seg_len_pool.append(curr_seg_len)
                    curr_seg_len = 0
            else:
                if curr_nucl == p[2].upper():
                    curr_len += 1
                    if p[2].upper() == p[2]:
                        curr_seg_len += 1
                    else:
                        seg_len_pool.append(curr_seg_len)
                        curr_seg_len=0
                else:
                    if curr_len > 2 and curr_len < 11:
                        hp_indel[curr_len - 3][curr_seg_len - curr_len] += 1
                            
                    curr_len = 1
                    curr_nucl = p[2].upper()
                    seg_len_pool = []
                    curr_seg_len = 0
                    if p[0] is not None:
                        i = p[0]
                        while i > -1 and read.seq[i] == curr_nucl:
                            curr_seg_len += 1
                            i -= 1
    
    print(f'Data: {dataset}')
                            
    for k in range(3, 11):
        x, y = [], []
        for key, value in hp_indel[k - 3].items():
            x.append(key)
            y.append(value)
        plt.plot(x, y)
        plt.title(f'Homopo len: {k}')
        plt.xlabel('difference with reference')
        plt.show()


# In[32]:



homo_quality_indel('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/b22_aln.bam', 'b22')


# In[33]:


homo_quality_indel('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/c24_aln.bam', 'c24')


# In[40]:


def errors_quality(filename, dataset):
    ins_quals = []
    subst_quals = []
    
    bamfile = pysam.AlignmentFile(filename, 'rb')
    for read in bamfile.fetch():
        idx_q = 0
        idx_aln = 0
        aln_pairs = read.get_aligned_pairs(with_seq=True)
        quals = read.query_qualities
        while idx_aln < len(aln_pairs):
            if aln_pairs[idx_aln][0] is None:
                idx_aln += 1
            elif aln_pairs[idx_aln][1] is None:
                ins_quals.append(quals[idx_q])
                idx_q += 1
                idx_aln += 1
            elif aln_pairs[idx_aln][2].lower() == aln_pairs[idx_aln][2]:
                subst_quals.append(quals[idx_q])
                idx_q += 1
                idx_aln += 1
            else: 
                idx_q += 1
                idx_aln += 1
                
    avg_ins_qual = sum(ins_quals) / len(ins_quals)
    avg_subst_qual = sum(subst_quals) / len(subst_quals)
    
    print('Data: {}\n'
          'Mean substitution: {:.3f}\n'
          'Mean insertion : {:.3f}\n'
          .format(dataset, avg_subst_qual,  avg_ins_qual))
    
    ins_distr_x, ins_distr_y = np.unique(ins_quals, return_counts=True)
    subst_distr_x, subst_distr_y = np.unique(subst_quals, return_counts=True)
    
    plt.figure(figsize=(7, 10))
    plt.subplot(211)
    plt.plot(ins_distr_x, ins_distr_y)
    plt.xlabel('quality')
    plt.title('insertion distribution')
    
    plt.subplot(212)
    plt.plot(subst_distr_x, subst_distr_y)
    plt.xlabel('quality')
    plt.title('substitution distribution')


# In[41]:


errors_quality('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/b22_aln.bam', 'b22')


# In[42]:


errors_quality('/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/c24_aln.bam', 'c24')


# In[48]:


filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/b22_bowtie_aln.bam'


print(showCover(filename))


# In[49]:


filename = '/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW2_2/Data/c24_bowtie_aln.bam'


print(showCover(filename))


# In[ ]:




