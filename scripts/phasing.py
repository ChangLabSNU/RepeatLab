import os
import sys
import re
import numpy as np
import pandas as pd
from itertools import chain
from sklearn.mixture import GaussianMixture
import scipy.stats as stats

SAMPLE_TARGET = sys.argv[1]
NUM_ALLELES = 2

sample_reads = {}

with open(f'/qbio/yoojung/repeat-expansion-diagnosis/analyses/{SAMPLE_TARGET}/result.log') as f:
    for line in f:
        if line.startswith('INFO: READ:'):
            fragments = line.split(' ')
            if fragments[2] == 'Status:True':
                repeat_count = int(fragments[3][11:])
                if repeat_count != 0:
                  read_id = fragments[1][5:]
                  sequence = fragments[6][7:]
                  pattern = fragments[5][9:]
                  sample_reads[read_id] = [repeat_count, sequence]

sample_reads_df = pd.DataFrame(sample_reads).T.reset_index().sort_values(by=[0])
sample_reads_df.columns = ['read_id', 'repeat_count', 'sequence']
sample_reads_df.reset_index(drop=True, inplace=True)

rc = sample_reads_df['repeat_count'].to_list()
sample_repeats = np.array(rc)
l10_sample_repeats = np.log10(sample_repeats)
sample_gmm = GaussianMixture(NUM_ALLELES).fit(l10_sample_repeats[:, np.newaxis])
gmm_labels = np.argsort(sample_gmm.means_.ravel())
gmmidx2label = {gmmidx: label for label, gmmidx in enumerate(gmm_labels)}
gmm_means = sample_gmm.means_[gmm_labels]
pred_labels = np.array([gmmidx2label[p] for p in sample_gmm.predict(l10_sample_repeats[:, np.newaxis])])
allele_repeats = 10 ** gmm_means

sample_reads_df.insert(2, 'allele', pred_labels)
sample_reads_df.insert(0, 'sample_id', SAMPLE_TARGET)

allele1_reads = sample_reads_df[sample_reads_df['allele'] == 0]['read_id'].to_list()
allele2_reads = sample_reads_df[sample_reads_df['allele'] == 1]['read_id'].to_list()

with open(f'/qbio/yoojung/repeat-expansion-diagnosis/on-target/phased_readID_list/{SAMPLE_TARGET}.allele1.readID.txt', 'w') as f:
    for read in allele1_reads:
        f.write(read + '\n')

with open(f'/qbio/yoojung/repeat-expansion-diagnosis/on-target/phased_readID_list/{SAMPLE_TARGET}.allele2.readID.txt', 'w') as f:
    for read in allele2_reads:
        f.write(read + '\n')