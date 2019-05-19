#!/usr/bin/env python
# coding: utf-8

# In[7]:


from Bio import Seq
import matplotlib.pyplot as plt
import time
import seaborn as sns

le = 100

u1 = []
u2 = []
u3 = []

for record in SeqIO.parse("/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW1/100x.1.fastq", "fastq"):

    c = record.seq
    d = record.letter_annotations["phred_quality"]

    u1.append((c.count('C') + c.count('G')) / 100)

    a2 = 0
    c2 = 0.001
    for n in range(0, le):
        if d[n] >= 20:
            c2 += 1
            if c[n] in {'G', 'C'}:
                a2 += 1
    u2.append(a2 / c2)

    if sum(i > 20 for i in d) > 80:
        u3.append((c.count('C') + c.count('G')) / 100)

plt.figure()
plt.title('frequency distribution % without any restrictions')
sns.distplot(u1)

plt.figure()
plt.title('frequency distribution % with quality>20 ')
sns.distplot(u2)
plt.figure()
sns.distplot(u3)


# In[8]:


from Bio import Seq
import matplotlib.pyplot as plt
import time
import seaborn as sns

le = 100

u1 = []
u2 = []
u3 = []

for record in SeqIO.parse("/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW1/100x.2.fastq", "fastq"):

    c = record.seq
    d = record.letter_annotations["phred_quality"]

    u1.append((c.count('C') + c.count('G')) / 100)

    a2 = 0
    c2 = 0.001
    for n in range(0, le):
        if d[n] >= 20:
            c2 += 1
            if c[n] in {'G', 'C'}:
                a2 += 1
    u2.append(a2 / c2)

    if sum(i > 20 for i in d) > 80:
        u3.append((c.count('C') + c.count('G')) / 100)

plt.figure()
plt.title('frequency distribution % without any restrictions')
sns.distplot(u1)

plt.figure()
plt.title('frequency distribution % with quality>20 ')
sns.distplot(u2)
plt.figure()
sns.distplot(u3)


# In[9]:


import numpy as np

def getQualityByPosition(file_name):
    posqv = []
    cnt = 0
    for rec in SeqIO.parse(file_name, "fastq"):
        if cnt % 100000 == 0:
            print(cnt)
        cnt += 1
        reclen = len(rec.seq)
        for i in range(reclen):
            if i >= len(posqv):
                posqv.append([])
            posqv[i].append(rec.letter_annotations["phred_quality"][i])
    return posqv

def showQualityByPosition(file_name):
    dataset = getQualityByPosition(file_name)
    fig = plt.figure(1, figsize=(30,30))
    ax = fig.add_subplot(111)
    ax.boxplot(dataset)
    plt.xlabel('position')
    plt.ylabel('quality')
    plt.title(file_name)
    plt.show()
    
showQualityByPosition("/Users/alena_paliakova/Google Drive/!Bioinf_drive/01_NGS/HW1/100x.1.fastq")


# In[ ]:




