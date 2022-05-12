import sys
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt

input_log = sys.argv[1]

repeat_count_pat = re.compile('\d+:\d+,')
sampleid_pat = re.compile('Input BAM file is .+.bam')
gene_pat = re.compile("p2sp .+',")
hist_dict = {}

with open(input_log) as logfile:
    for line in logfile:
        repeat_count = re.findall(repeat_count_pat, line)
        if len(repeat_count) != 0:
            count_li = [c.strip(',') for c in repeat_count]
            for s in count_li:
                count_set = s.split(':')
                repeat_len = int(count_set[0])
                read_num = int(count_set[1])
                hist_dict[repeat_len] = read_num
                
        sampleid_sen = re.findall(sampleid_pat, line)
        if len(sampleid_sen) != 0:
            sampleid = str(sampleid_sen[0]).split('/')[-1].split('_')[0]
            
        gene_sen = re.findall(gene_pat, line)
        if len(gene_sen) != 0:
            gene = gene_sen[0].split("'")[1].upper()
            
fig, axes = plt.subplots(1,2,figsize=(12,7))

ax1 = axes[0]
# ax1.bar(hist_dict.keys(), hist_dict.values(), color='#22577E')
ax1.hist(hist_dict, color='#22577E')
ax1.set_xlabel('repeat count\n(length)')
ax1.set_ylabel('number of reads\n(proportion)')
ax1.set_title('{}_{}'.format(sampleid, gene))

ax2 = axes[1]
column_labels = ['repeat count', 'number of reads']
arraydata = np.array([i for i in hist_dict.items()])
width = [0.4]*(len(arraydata)+1)
ax2.axis('off')
ax2.table(cellText=arraydata, colLabels=column_labels, colColours=['#E4E9BE']*3, loc='center', cellLoc='center', colWidths=width).set_fontsize(8)

fig.savefig('{}_{}_RepeatCount.png'.format(sampleid, gene))