#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 10:58:04 2019

@author: annie
"""

import numpy as np
import gzip
import matplotlib.pyplot as plt

files = ['/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz','/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz','/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz','/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz']

# files = ['I_HATE_TEST_FILE.txt.gz','test2.txt.gz']

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter)-33

for k in range(len(files)): 
    with gzip.open(files[k],'rt') as fh:
        line = fh.readline()
        line = fh.readline()
        line = line.strip('\n')
        # print(line)
        bp = len(line)
        scores = np.zeros(bp)
        i = 2
        for line in fh:  
            i +=1
            line = line.strip('\n')
            # if i ==2:
            #     bp = len(line)
            #     scores = np.zeros(bp)

            if i % 4 ==0:
                # print(line)
                for j in range(bp):
                    scores[j]+=convert_phred(line[j])
            if i % 50000000 ==0:
                print(i)
        meanscores = scores*4/i

        plt.bar(range(1,bp+1),meanscores)
        plt.xlabel('Base position')
        plt.ylabel('Mean quality scores')
        plt.savefig('mqdis_1294_S1_L008_R'+str(k+1)+'_001'+'.png')
        plt.close()
