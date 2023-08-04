import sys
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

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

def GMM_clustering(l10_sample_repeats, NUM_ALLELES):
    sample_gmm = GaussianMixture(NUM_ALLELES, random_state=1).fit(l10_sample_repeats[:, np.newaxis])
    gmm_labels = np.argsort(sample_gmm.means_.ravel())
    gmmidx2label = {gmmidx: label for label, gmmidx in enumerate(gmm_labels)}
    gmm_means = sample_gmm.means_[gmm_labels]
    pred_labels = np.array([gmmidx2label[p] for p in sample_gmm.predict(l10_sample_repeats[:, np.newaxis])])

    return pred_labels

pred_labels = GMM_clustering(l10_sample_repeats, NUM_ALLELES)
label_matched_repeats = []
for repeat, label in zip(l10_sample_repeats, pred_labels):
    label_matched_repeats.append([repeat, label])
round_num = 1

while list(pred_labels).count(0) / len(pred_labels) < 0.1 or list(pred_labels).count(1) / len(pred_labels) < 0.1:
    if list(pred_labels).count(0) / len(pred_labels) < 0.1:
        l10_sample_repeats = np.array([repeat for repeat, label in label_matched_repeats if label == 1])
    elif list(pred_labels).count(1) / len(pred_labels) < 0.1:
        l10_sample_repeats = np.array([repeat for repeat, label in label_matched_repeats if label == 0])
    pred_labels = GMM_clustering(l10_sample_repeats, NUM_ALLELES)
    label_matched_repeats = []
    for repeat, label in zip(l10_sample_repeats, pred_labels):
        label_matched_repeats.append([repeat, label])
    round_num += 1
    if round_num > 3:
        break

if len(set(pred_labels)) == 1:
    pred_labels = np.array([0] * (len(pred_labels)//2) + [1] * (len(pred_labels) - len(pred_labels)//2))

removed = [l for l in rc if l not in list(np.around(10**l10_sample_repeats).astype(int))]
sample_reads_df = sample_reads_df[sample_reads_df['repeat_count'].isin(removed) == False]

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