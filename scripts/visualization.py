import os
import sys
import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.io as pio
import re 
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import statsmodels.api as sm
import mpld3
from mpld3 import plugins
from sklearn.mixture import GaussianMixture
from scipy import optimize
from itertools import chain
from Bio import SeqIO
from collections import defaultdict

ANALYSIS = sys.argv[1]
SAMPLE_NAME = ANALYSIS.split('-')[0]
TARGET_GENE = ANALYSIS.split('-')[1]
OUTPUT_DIR_T = f'analyses/{ANALYSIS}'

seq_summary = None

def process_form_summary(file_path):
    seq_summary = pd.read_csv(
        file_path, delimiter='\t',
        usecols=[
            'channel', 'start_time', 'duration',
            'sequence_length_template', 'mean_qscore_template'])
    seq_summary.sort_values('start_time', inplace=True)
    seq_summary.reset_index(drop=True, inplace=True)
    seq_summary = seq_summary
    return seq_summary

rough_basecall_summary = f'basecalls/dorado/fast/{SAMPLE_NAME}/{SAMPLE_NAME}.fast-called.summary.txt'
seq_summary = process_form_summary(rough_basecall_summary)

pattern_file = 'pre-requisites/RepeatHMM/bin/reference_sts/hg38/hg38.predefined.pa'
target_gene_lower = TARGET_GENE.lower()

with open('pre-requisites/RepeatHMM/bin/reference_sts/hg38/hg38.predefined.pa') as f:
    for line in f:
        if target_gene_lower in line:
            target_region_info = line.split(',')[:-1]
            break

if not os.path.exists(f'{OUTPUT_DIR_T}/others'):
    os.mkdir(f'{OUTPUT_DIR_T}/others')

with open(f'{OUTPUT_DIR_T}/others/target_region.bed', 'w') as f:
    target_region_df = pd.DataFrame(target_region_info, index=['gene', 'chr', 'start', 'end', 'pattern', 'strand&range', 'others']).T
    gene = target_region_df.loc[0, 'gene']
    chr = target_region_df.loc[0, 'chr']
    start = target_region_df.loc[0, 'start']
    end = target_region_df.loc[0, 'end']
    f.write('{}\t{}\t{}\t{}'.format(chr, start, end, gene))

region_bed = None
target_names = None

def process_form_bed(file_path):

    region_bed = pd.read_csv(
        file_path, sep='\t', header=None,
        names=[
            'chrom', 'start', 'end', 'tname'])

    region_bed['region'] = [
        '{}:{}-{}'.format(x['chrom'], x.start, x.end)
        for _, x in region_bed.iterrows()]
    target_names = region_bed['tname'].to_list()

    return region_bed, target_names

region_bed, target_names = process_form_bed(f'{OUTPUT_DIR_T}/others/target_region.bed')

# Basecall QC - Pass / Fail reads proportion plot

seq_summary['passes_filtering'] = np.where(seq_summary['mean_qscore_template'] >= 6, True, False)
seq_summary.passes_filtering = seq_summary.passes_filtering.astype('category')
pass_fail = seq_summary.groupby('passes_filtering', observed=False).size()

pass_percentage = 100 * pass_fail[True] / len(seq_summary)
fail_percentage = 100 - pass_percentage

basecall_QC_PF_fig = go.Figure(data=[
    go.Bar(
        y=['Total'],
        x=[pass_percentage],
        orientation='h',
        name='Pass',
        marker=dict(color='#54B8B1')
    ),
    go.Bar(
        y=['Total'],
        x=[fail_percentage],
        orientation='h',
        name='Fail',
        marker=dict(color='#EF4135')
    )
])

basecall_QC_PF_fig.update_layout(
    title="Basecalling : Pass / Fail reads",
    xaxis=dict(title='%age Reads', range=[0, 100]),
    barmode='stack',
    plot_bgcolor="#f4f4f4"
)

basecall_QC_PF = pio.to_html(basecall_QC_PF_fig, full_html=False)

# Read length distribution plot

sorted_lengths = seq_summary.sequence_length_template.sort_values(ascending=False).reset_index(drop=True)
cumulative_length = sorted_lengths.cumsum()
total_bases = cumulative_length.iloc[-1]
mean_length = total_bases / len(seq_summary)
n50_index = cumulative_length.searchsorted(total_bases / 2)
n50_length = sorted_lengths.iloc[n50_index]

LOGPLOT_FOCUS_PCT = 2
LOGPLOT_FLANKING_VIEW = 0.05

datas = [seq_summary.sequence_length_template]

ax_focus = np.percentile(datas[0], [LOGPLOT_FOCUS_PCT, 100 - LOGPLOT_FOCUS_PCT])
ax_focuswidth = ax_focus[1] - ax_focus[0]
ax_view = (
    ax_focus[0] - ax_focuswidth * LOGPLOT_FLANKING_VIEW,
    ax_focus[1] + ax_focuswidth * LOGPLOT_FLANKING_VIEW
)

basecall_QC_readDist_fig = go.Figure()

basecall_QC_readDist_fig.add_trace(go.Histogram(x=datas[0], nbinsx=400, marker_color='#0098A9', marker_line_color='black', marker_line_width=1))

basecall_QC_readDist_fig.update_layout(
    title="Read length distribution",
    xaxis=dict(title='Read Length (bases)', range=[ax_view[0], ax_view[1]]),
    yaxis=dict(title='Number of reads'),
    shapes=[
        {
            'type': 'line',
            'x0': mean_length,
            'x1': mean_length,
            'y0': 0,
            'y1': 1,
            'line': {'color': 'black', 'width': 1, 'dash': 'dash'},
        },
        {
            'type': 'line',
            'x0': n50_length,
            'x1': n50_length,
            'y0': 0,
            'y1': 1,
            'line': {'color': 'black', 'width': 1, 'dash': 'dash'},
        }
    ],
    annotations=[
        {
            'x': mean_length,
            'y': 1,
            'xref': 'x',
            'yref': 'paper',
            'text': 'Mean: {:.0f}'.format(mean_length),
            'showarrow': True,
            'arrowhead': 2,
            'arrowwidth': 2,
            'arrowcolor': '#636363',
            'ax': 0,
            'ay': -40,
        },
        {
            'x': n50_length,
            'y': 1,
            'xref': 'x',
            'yref': 'paper',
            'text': 'N50: {}'.format(n50_length),
            'showarrow': True,
            'arrowhead': 2,
            'arrowwidth': 2,
            'arrowcolor': '#636363',
            'ax': 0,
            'ay': -40,
        }
    ],
    plot_bgcolor="#f4f4f4"
)

basecall_QC_readDist = pio.to_html(basecall_QC_readDist_fig, full_html=False)

# Repeat count plot
NUM_ALLELES = 2

input_log = f'{OUTPUT_DIR_T}/result.log'
p2hmm_pat = re.compile('INFO: *p2bamhmm (.*)')
logcontent = open(input_log).read()

# Get the sample id
sampleid = SAMPLE_NAME

# Get the gene name and repeat counts
p2hmm_found = p2hmm_pat.search(logcontent)
p2hmminfo = eval(p2hmm_found.groups()[0])
gene = TARGET_GENE
repcounts = [
    tuple(map(int, inst.rstrip(',').split(':')))
    for inst in p2hmminfo[3].split(':', 1)[1].split(' ')]
rc = list(chain(*[[length] * count for length, count in repcounts]))
rc = [r for r in rc if r != 0]   # remove 0

# Fit a Gaussian mixture model to the loaded repeat counts
sample_repeats = np.array(rc)
l10_sample_repeats = np.log10(sample_repeats)

def GMM_clustering(l10_sample_repeats, NUM_ALLELES):
    sample_gmm = GaussianMixture(NUM_ALLELES, random_state=1).fit(l10_sample_repeats[:, np.newaxis])
    gmm_labels = np.argsort(sample_gmm.means_.ravel())
    gmmidx2label = {gmmidx: label for label, gmmidx in enumerate(gmm_labels)}
    gmm_means = sample_gmm.means_[gmm_labels]
    pred_labels = np.array([gmmidx2label[p] for p in sample_gmm.predict(l10_sample_repeats[:, np.newaxis])])

    return sample_gmm, gmm_labels, gmmidx2label, gmm_means, pred_labels

sample_gmm, gmm_labels, gmmidx2label, gmm_means, pred_labels = GMM_clustering(l10_sample_repeats, NUM_ALLELES)

label_matched_repeats = []
for repeat, label in zip(l10_sample_repeats, pred_labels):
    label_matched_repeats.append([repeat, label])
round_num = 1

while list(pred_labels).count(0) / len(pred_labels) < 0.1 or list(pred_labels).count(1) / len(pred_labels) < 0.1:
  if list(pred_labels).count(0) / len(pred_labels) < 0.1:
    l10_sample_repeats = np.array([repeat for repeat, label in label_matched_repeats if label == 1])
  elif list(pred_labels).count(1) / len(pred_labels) < 0.1:
    l10_sample_repeats = np.array([repeat for repeat, label in label_matched_repeats if label == 0])
  sample_gmm, gmm_labels, gmmidx2label, gmm_means, pred_labels = GMM_clustering(l10_sample_repeats, NUM_ALLELES)
  label_matched_repeats = []
  for repeat, label in zip(l10_sample_repeats, pred_labels):
    label_matched_repeats.append([repeat, label])
  round_num += 1
  if round_num > 3:
    break

if len(set(pred_labels)) == 1:
    pred_labels = np.array([0] * (len(pred_labels)//2) + [1] * (len(pred_labels) - len(pred_labels)//2))

removed = [l for l in rc if l not in list(np.around(10**l10_sample_repeats).astype(int))]
sample_repeats = list(np.around(10**l10_sample_repeats).astype(int))

# Plotting
BINS = 25
ALLELE_SHADE_COLORS = ['#FBC5C5', '#C7D36F']
ALLELE_MARKER_COLORS = ['#ff0000', '#008000']

LOGPLOT_FOCUS_PCT = 2
LOGPLOT_BAR_MARGIN = 0.1
LOGPLOT_FLANKING_VIEW = 0.05
LOGPLOT_DROP_SPINES = 5
LOGPLOT_PDF_ADDITIONAL_SCALE = 1

fig, ax = plt.subplots(1, 1, figsize=(7, 3))

plt.suptitle('{}-{}'.format(sampleid, gene), y=1.3, fontsize=23, fontweight='bold')
plt.subplots_adjust(wspace=0.4)

plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor']='white'

ax.set_xscale('log')

# Set up the ranges to compute histogram and to show
ax_focus = np.percentile(l10_sample_repeats, [LOGPLOT_FOCUS_PCT, 100 - LOGPLOT_FOCUS_PCT])
ax_focuswidth = ax_focus[1] - ax_focus[0]
ax_view = (
    ax_focus[0] - ax_focuswidth * LOGPLOT_FLANKING_VIEW,
    ax_focus[1] + ax_focuswidth * LOGPLOT_FLANKING_VIEW)

# Plot the histogram for each allele and annotate with the stats
prevcount = np.zeros(BINS)
totalreads = len(sample_repeats)
allele_labels = []
for i, shadecolor, markercolor in zip(range(NUM_ALLELES), ALLELE_SHADE_COLORS, ALLELE_MARKER_COLORS):
    l10counts_allele = l10_sample_repeats[pred_labels == i]
    mediancount = np.around(10 ** np.median(l10counts_allele))

    # Plot histogram
    count, left = np.histogram(l10counts_allele, bins=BINS, range=ax_focus)
    margin = (left[1] - left[0]) * LOGPLOT_BAR_MARGIN
    barleft = 10 ** (left[:-1] + margin/2)
    ax.bar(barleft, count,
            width=10 ** (left[1:] - margin) - barleft,
            color=shadecolor, bottom=prevcount, edgecolor=markercolor, linewidth=.5)
    prevcount = count + prevcount
    
    def add_tooltip(ax, x, y, labels):
        points = ax.plot(x, y, 'o', markersize=5)

        tooltip = plugins.PointHTMLTooltip(points[0], labels)
        plugins.connect(fig, tooltip)

    labels = [f'Allele #{i+1}<br>median={mediancount:g}<br>{len(l10counts_allele):} reads']
    allele_labels.append(f'Allele #{i+1} : repeat count = {mediancount:g}, total {len(l10counts_allele):} reads')
    add_tooltip(ax, mediancount, np.max(prevcount)/2, labels)

# Find the optimal scaling factor for the dist plot to overlay under the histogram
barcenters = (left[1:] + left[:-1]) / 2
barprobs = np.exp(sample_gmm.score_samples(barcenters[:, None]))
def evaluate_scaling(x):
    return np.sum((prevcount - x * barprobs) ** 2)
prob_scale = optimize.minimize(evaluate_scaling, (1,)).x[0]

# Plot the probability distribution of the GMM
pdf_x = np.linspace(ax_view[0], ax_view[1], 200)
pdf_y = sample_gmm.score_samples(pdf_x.reshape(-1, 1))
ax.plot(10 ** pdf_x, np.exp(pdf_y) * prob_scale * LOGPLOT_PDF_ADDITIONAL_SCALE,
         color='#432686', linewidth=1.5)

# Adjust the aesthetics
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_position(('outward', LOGPLOT_DROP_SPINES))
ax.spines['bottom'].set_position(('outward', LOGPLOT_DROP_SPINES))

ax.grid(which='major', alpha=0.5)
ax.grid(which='minor', alpha=0.2)

ax.set_xlabel('Repeat count', fontsize=18, labelpad=15)
ax.set_ylabel('Read count', fontsize=18, labelpad=15)

# Show the tick labels in the regular number format (not exponential)
if 10**ax_view[1] - 10**ax_view[0] < 100:
  ax.set_xticks(np.arange(np.floor(10**ax_view[0])+1, np.ceil(10**ax_view[1]),np.ceil(np.around(10**ax_view[1] - 10**ax_view[0])/10)+1))
  ticks = ax.get_xticks()
  ax.set_xticklabels(list(map(int, ticks)))
else:
  ticks = ax.get_xticks()
  ax.set_xticks(ticks)
  ax.set_xticklabels(list(map(int, ticks)))

repeat_count_plot = mpld3.fig_to_html(fig)

if len(removed) > 0:
  print(f'** NOTICE (outlier existence)\n : {len(removed)} read(s) of {", ".join([str(num) for num in removed])} repeat counts is(are) considered as carryover contamination and eliminated.')

# Format the raw count data for display
output_df = pd.Series(sample_repeats).value_counts().sort_index().reset_index()
output_df.columns = ['Repeat size', 'Read count']
output_df['Allele'] = [
    f'#{gmmidx2label[l] + 1}'
    for l in sample_gmm.predict(np.log10(output_df['Repeat size'].values.reshape(-1, 1)))]
output_df.reset_index(drop=True, inplace=True)
output_df[['Allele', 'Repeat size', 'Read count']]

# Repeat sequence visualization
sample_reads = {}

with open(f'{OUTPUT_DIR_T}/result.log') as f:
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

with open(f'{OUTPUT_DIR_T}/others/repeat_seq.fasta', 'w') as f:
    for i in range(len(sample_reads_df)):
        f.write('>{}\n{}\n'.format('Read_'+str(i+1), sample_reads_df.iloc[i]['sequence']))

REPEAT_COLORS = ['bright-cyan', 'bright-red', 'bright-green', 'bright-yellow', 'bright-blue', 'bright-magenta']

def profile_tandem_repeats(seq, kmersize=2):
  kmercounts = defaultdict(int)
  for i in range(len(seq) - kmersize + 1):
    kmer = seq[i:i+kmersize]
    if '^' not in kmer: # Skip over repeat marks
      kmercounts[kmer] += 1
  return kmercounts

pat_repmarks = re.compile('\\([^)]*\\)')

def summarize_tandem_repeats(seq, kmer_short=3, kmer_long=6, minscore=6):
  candidates = []

  seq_wo_repmarks = pat_repmarks.sub('^', seq)
  for ksize in range(kmer_short, kmer_long + 1):
    kmercounts = profile_tandem_repeats(seq_wo_repmarks, ksize)
    for kmer, count in kmercounts.items():
      if count >= minscore:
        repoccurrences = re.findall(f'(({kmer}){{2,}})', seq_wo_repmarks)
        totalreplength = sum(len(full) for full, repunit in repoccurrences)
        if totalreplength >= minscore:
          candidates.append([totalreplength, ksize, kmer])

  if not candidates:
    return seq

  candidates.sort(key=lambda x: (-x[0], x[1]))

  # Replace the occurrences of the first candidate with repeat marks
  repeat = candidates[0][2]
  def replace_repeat(match):
    repcount = len(match.group(0)) // len(repeat)
    if repcount <= 1:
      return match.group(0)
    else:
      return f'({repeat}/{repcount})'

  replaced = re.sub(f'({repeat})+', replace_repeat, seq)
  if replaced == seq:
    return seq
  else:
    return summarize_tandem_repeats(replaced, kmer_short, kmer_long, minscore)

pat_repmarks_full = re.compile('(\\(([^)]*)/(\\d+)\\)|[^()/]+)')

def colorize_repeats(seq):
  # Find all repeats and their occurrences
  repcounts = defaultdict(int)
  for m in pat_repmarks_full.finditer(seq):
    rep, repunit, count = m.groups()
    if repunit is not None:
      repcounts[repunit] += int(count)

  repeats = sorted(repcounts.items(), key=lambda x: (-x[1], x[0]))
  if len(repeats) > len(REPEAT_COLORS):
    REPEAT_COLORS[:] = REPEAT_COLORS * 2 # reuse color scheme for another cycle
  repeatcolors = {rep: color for (rep, cnt), color in zip(repeats, REPEAT_COLORS)}

  outputs = []
  for m in pat_repmarks_full.finditer(seq):
    rep, repunit, count = m.groups()
    if repunit is None or repunit not in repeatcolors:
      outputs.append('<span class="gray">' + rep + '</span>')
    else:
        if int(count) >= 10:
            color = repeatcolors[repunit]
            outputs.append(f'<span class="{color}">({repunit}/{count})</span>')
        else:
            color = repeatcolors[repunit]
            outputs.append(f'<span class="{color}">({repunit}/{count})</span>')
  return ''.join(outputs)

turnover = np.bincount(pred_labels)[0]

if len(rc) > 30:
  subsampled_readnum = np.random.choice(len(rc), 30, replace=False)
else:
  subsampled_readnum = np.array([i+1 for i in range(len(rc)-1)])

rs_line = '-'*50
rs_allele1_label = '|' + ' '*19 + ' Allele 1 ' + ' '*19+ '|'
rs_allele2_label = '|' + ' '*19 + ' Allele 2 ' + ' '*19 + '|'

# Methylation plot

allele1_meth_prof = pd.read_csv(f'meth-profiles/{ANALYSIS}.allele1.meth-called.bedMethyl', sep='\t', names=['chr', 'start_pos', 'end_pos', 'modified_base_code', 'score', 'strand', 'start', 'end', 'color', 'N_numbers'])
allele1_N_numbers = allele1_meth_prof['N_numbers'].str.split(' ', expand=True)
allele1_N_numbers.columns = ['coverage', 'fraction_modified', 'N_mod', 'N_canonical', 'N_other_mod', 'N_delete', 'N_fail', 'N_diff', 'N_nocall']
allele1_meth_prof = pd.concat([allele1_meth_prof, allele1_N_numbers], axis=1)
allele1_meth_prof = allele1_meth_prof.astype({'coverage' : 'int', 'fraction_modified' : 'float'})

allele2_meth_prof = pd.read_csv(f'meth-profiles/{ANALYSIS}.allele2.meth-called.bedMethyl', sep='\t', names=['chr', 'start_pos', 'end_pos', 'modified_base_code', 'score', 'strand', 'start', 'end', 'color', 'N_numbers'])
allele2_N_numbers = allele2_meth_prof['N_numbers'].str.split(' ', expand=True)
allele2_N_numbers.columns = ['coverage', 'fraction_modified', 'N_mod', 'N_canonical', 'N_other_mod', 'N_delete', 'N_fail', 'N_diff', 'N_nocall']
allele2_meth_prof = pd.concat([allele2_meth_prof, allele2_N_numbers], axis=1)
allele2_meth_prof = allele2_meth_prof.astype({'coverage' : 'int', 'fraction_modified' : 'float'})

allele1_meth_prof, allele2_meth_prof = allele1_meth_prof[allele1_meth_prof['coverage'] > 5], allele2_meth_prof[allele2_meth_prof['coverage'] > 5]

DISCRETE_COLOR_LEVELS = 10
_cmap = mpl.colormaps['YlGnBu']
_cmap.set_gamma(0.6)
sitemarkers_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'mycolormap', _cmap(np.linspace(0, 1, DISCRETE_COLOR_LEVELS)),
    DISCRETE_COLOR_LEVELS)

ALLELE1_COLOR = '#444444'
ALLELE2_COLOR = '#e64980'

def plot_methylation_with_trend(ax, meth_prof, color, marker, allele, repeat_region):
  ax.scatter(meth_prof['start_pos'], meth_prof['fraction_modified'], zorder=6,
              edgecolor='none', facecolor=color,
              s={'s': 15, 'o': 20}[marker], alpha=.8,
              marker=marker, label=allele)

  trend = sm.nonparametric.lowess(exog=meth_prof['start_pos'], endog=meth_prof['fraction_modified'], frac=0.4)
  ax.plot(trend[:,0], trend[:,1], c=color, linewidth=1.5, linestyle='-', zorder=4)
  ax.set_xticks(ax.get_xticks())
  ax.set_yticks(np.arange(0,120,20))
  ax.set_xticklabels((ax.get_xticks()/1000).astype(int), size=12)
  ax.set_yticklabels(np.arange(0,120,20), size=12)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_position(('outward', LOGPLOT_DROP_SPINES))
  ax.spines['bottom'].set_position(('outward', LOGPLOT_DROP_SPINES))
  ax.legend()
  target_chr = repeat_region.split(':')[0]
  repeat_window = repeat_region.split(':')[1].split('-')
  repeat_start, repeat_end = int(repeat_window[0]), int(repeat_window[1])
  repeat_length = repeat_end - repeat_start
  ax.set_xlabel(f'Genomic coordinate in {target_chr} (x1000)', size=12)
  ax.set_ylabel('Methylation rate (%)', size=12)

meth_fig, ax = plt.subplots(1,1, figsize=(18,3))

target_region = str(region_bed.loc[0, 'region'])
plot_methylation_with_trend(ax, allele1_meth_prof, ALLELE1_COLOR, 's', 'allele 1', target_region)
plot_methylation_with_trend(ax, allele2_meth_prof, ALLELE2_COLOR, 'o', 'allele 2', target_region)

methylation_plot = mpld3.fig_to_html(meth_fig)

a1_keys_meth_mean, a2_keys_meth_mean = None, None
meth_key_df = None

def plot_methylation_key_region(ax, allele1_meth_prof, allele2_meth_prof, target_gene):

  if target_gene == 'DMPK':

      methylation_key_positions = "45768652,45768667,45768673,45768678,45768682,45768687 45770725,45770739,45770744,45770750,45770784,45770788 45769906,45769912,45769924,45769933,45769951,45769953"
      CTCF1_positions = "45770307,45770328,45770332,45770342,45770348,45770371,45770381,45770385,45770390,45770408,45770415,45770417,45770429,45770436,45770455,45770463,45770469,45770495,45770497,45770501,45770512,45770515,45770525,45770537,45770540"
      CTCF2_positions = "45769994,45770013,45770015,45770022,45770039,45770068,45770076,45770091,45770107,45770111,45770116"

      meth_key_groups = methylation_key_positions.split(" ")
      meth_keys = []
      for group in meth_key_groups:
          meth_keys.append(list(map(int, group.split(","))))

      CTCF1_sites = [int(site) for site in CTCF1_positions.split(",")]
      CTCF2_sites = [int(site) for site in CTCF2_positions.split(",")]

      a1_CTCF1_meth_mean = np.mean(allele1_meth_prof.loc[allele1_meth_prof['start_pos'].isin(CTCF1_sites), 'fraction_modified'].tolist())
      a1_CTCF2_meth_mean = np.mean(allele1_meth_prof.loc[allele1_meth_prof['start_pos'].isin(CTCF2_sites), 'fraction_modified'].tolist())
      a1_g1_meth_mean = np.mean(allele1_meth_prof.loc[allele1_meth_prof['start_pos'].isin(meth_keys[0]), 'fraction_modified'].tolist())
      a1_g2_meth_mean = np.mean(allele1_meth_prof.loc[allele1_meth_prof['start_pos'].isin(meth_keys[1]), 'fraction_modified'].tolist())
      a1_g3_meth_mean = np.mean(allele1_meth_prof.loc[allele1_meth_prof['start_pos'].isin(meth_keys[2]), 'fraction_modified'].tolist())

      a2_CTCF1_meth_mean = np.mean(allele2_meth_prof.loc[allele2_meth_prof['start_pos'].isin(CTCF1_sites), 'fraction_modified'].tolist())
      a2_CTCF2_meth_mean = np.mean(allele2_meth_prof.loc[allele2_meth_prof['start_pos'].isin(CTCF2_sites), 'fraction_modified'].tolist())
      a2_g1_meth_mean = np.mean(allele2_meth_prof.loc[allele2_meth_prof['start_pos'].isin(meth_keys[0]), 'fraction_modified'].tolist())
      a2_g2_meth_mean = np.mean(allele2_meth_prof.loc[allele2_meth_prof['start_pos'].isin(meth_keys[1]), 'fraction_modified'].tolist())
      a2_g3_meth_mean = np.mean(allele2_meth_prof.loc[allele2_meth_prof['start_pos'].isin(meth_keys[2]), 'fraction_modified'].tolist())

      a1_keys_meth_mean = np.array([a1_CTCF1_meth_mean, a1_CTCF2_meth_mean, a1_g1_meth_mean, a1_g2_meth_mean, a1_g3_meth_mean])
      a2_keys_meth_mean = np.array([a2_CTCF1_meth_mean, a2_CTCF2_meth_mean, a2_g1_meth_mean, a2_g2_meth_mean, a2_g3_meth_mean])

  elif os.path.exists(f'key_sites.txt'):
      with open(f'key_sites.txt', 'r') as f:
          if len(f.read().splitlines()) > 0:
              key_sites = f.read().splitlines()[0]
              meth_key_groups = key_sites.split(" ")
              if len(meth_key_groups) > 0:
                  meth_keys = []
                  a1_keys_meth_mean = []
                  a2_keys_meth_mean = []
                  for group in meth_key_groups:
                      meth_keys.append(list(map(int, group.split(","))))
                  for key in meth_keys:
                      a1_keys_meth_mean.append(np.mean(allele1_meth_prof.loc[allele1_meth_prof['start_pos'].isin(key), 'fraction_modified'].tolist()))
                      a2_keys_meth_mean.append(np.mean(allele2_meth_prof.loc[allele2_meth_prof['start_pos'].isin(key), 'fraction_modified'].tolist()))
                  a1_keys_meth_mean = np.array(a1_keys_meth_mean)
                  a2_keys_meth_mean = np.array(a2_keys_meth_mean)
          else:
              return print('No methylation key regions.')

  else:
      return print('No methylation key regions.')

  if len(a1_keys_meth_mean) > 0 and len(a2_keys_meth_mean) > 0:
    from matplotlib.patches import Polygon
    from matplotlib.colors import Normalize
    import matplotlib.cm as cm

    DISCRETE_COLOR_LEVELS = 10
    _cmap = plt.get_cmap('YlGnBu').copy()
    _cmap.set_gamma(0.6)
    _cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'mycolormap', _cmap(np.linspace(0, 1, DISCRETE_COLOR_LEVELS)),
        DISCRETE_COLOR_LEVELS)
    cmap_top = _cmap
    cmap_bottom = _cmap

    norm_top = Normalize(vmin=0, vmax=100)
    norm_bottom = Normalize(vmin=0, vmax=100)

    gap = 0.05

    for i in range(len(a1_keys_meth_mean)):
        x_left, x_right = i+gap, i+1-gap
        y_bottom, y_top = 0, 1

        coords_top_left = [(x_left, y_top), (x_left, y_bottom), (x_right, y_bottom)]
        coords_bottom_right = [(x_right, y_top), (x_left, y_top), (x_right, y_bottom)]

        poly_top = Polygon(coords_top_left, closed=True)
        poly_bottom = Polygon(coords_bottom_right, closed=True)

        val_top = a1_keys_meth_mean[i]
        val_bottom = a2_keys_meth_mean[i]

        poly_top.set_facecolor(cmap_top(norm_top(val_top)))
        poly_bottom.set_facecolor(cmap_bottom(norm_bottom(val_bottom)))

        poly_top.set_edgecolor('none')
        poly_bottom.set_edgecolor('none')

        ax.add_patch(poly_top)
        ax.add_patch(poly_bottom)

        linewidth=1
        ax.plot([i+gap, i+1-gap], [1, 0], color='black', linewidth=linewidth, zorder=3)

        ax.plot([i+gap, i+gap], [0, 1], color='black', linewidth=linewidth, zorder=3)
        ax.plot([i+1-gap, i+1-gap], [0, 1], color='black', linewidth=linewidth, zorder=3)
        ax.plot([i+gap, i+1-gap], [0, 0], color='black', linewidth=linewidth, zorder=3)
        ax.plot([i+gap, i+1-gap], [1, 1], color='black', linewidth=linewidth, zorder=3)

    ax.set_xlim(0, len(a1_keys_meth_mean))
    ax.invert_yaxis()
    ax.set_yticks([])

    for axis in ['top', 'bottom', 'left', 'right']:
      ax.spines[axis].set_visible(False)

    sm_top = cm.ScalarMappable(norm=norm_top, cmap=cmap_top)
    cbar = plt.colorbar(sm_top, ax=ax, location='bottom', pad=0.3, aspect=15, shrink=0.6)
    cbar.set_label('Methylation rate (%)', fontsize=10)

    ax.set_xticks(np.arange(len(a1_keys_meth_mean)) + 0.5)
    if target_gene == 'DMPK':
      labels = ['CTCF1', 'CTCF2', 'Locus1', 'Locus2', 'Locus3']
    else:
      labels = [f'Locus{i+1}' for i in range(len(a1_keys_meth_mean))]
    ax.set_xticklabels(labels, fontsize=15)
    meth_key_df = pd.DataFrame([a1_keys_meth_mean, a2_keys_meth_mean], columns=labels, index=['Allele1', 'Allele2'])
    return meth_key_df

meth_keys_fig, ax = plt.subplots(1, 1, figsize=(6, 2))

meth_key_df = plot_methylation_key_region(ax, allele1_meth_prof, allele2_meth_prof, TARGET_GENE)
if TARGET_GENE == 'DMPK':
  labels = ['CTCF1', 'CTCF2', 'Locus1', 'Locus2', 'Locus3']
else:
  labels = [f'Locus{i+1}' for i in range(len(a1_keys_meth_mean))]
methylation_key_plot = mpld3.fig_to_html(meth_keys_fig)

with open(f'{OUTPUT_DIR_T}/result.html', 'w') as f:
    f.write('<!DOCTYPE html>')
    f.write('<html lang="en">')
    f.write('<head>')
    f.write('<meta charset="UTF-8">')
    f.write('<meta name="viewport" content="width=device-width, initial-scale=1.0">')
    f.write('<title>Data Report ({})</title>'.format(ANALYSIS))
    f.write('<link href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css" rel="stylesheet">')
    f.write('<style> .gray { color: #808080; } .bright-cyan { color: #36428a; } .bright-red { color: #ff0000; } .bright-green { color: #1d591f; } .bright-yellow { color: #ffff00; } .bright-blue { color: #0000ff; } .bright-magenta { color: #ff00ff; }</style>')
    f.write('</head>')
    f.write('<body>')
    f.write('<div class="container">')
    f.write('<h1>Data Report ({})</h1>'.format(ANALYSIS))
    f.write('<h1>1. Basecall QC</h1>')
    f.write('<h2>Basecall Pass / Fail reads proportion</h2>')
    f.write('<h6>The total number of reads is: {}</h6>'.format(len(seq_summary)))
    f.write('<h6>The number of basecall passed reads is: {}</h6>'.format(pass_fail[True]))
    f.write(f'{basecall_QC_PF}')
    f.write('<h2>Read Length Distribution</h2>')
    f.write(f'{basecall_QC_readDist}')
    f.write('<h1>2. Repeat Count</h1>')
    f.write('<h2>Repeat count plot</h2>')
    f.write('<h6>{}</h6>'.format(allele_labels[0]))
    f.write('<h6>{}</h6>'.format(allele_labels[1]))
    f.write(repeat_count_plot)
    f.write('<h2>Repeat count table</h2>')
    f.write(output_df.to_html())
    f.write('<h1>3. Repeat Sequence</h1>')
    f.write('<h6>{}</h6>'.format(rs_line))
    f.write('<h6>{}</h6>'.format(rs_allele1_label))
    f.write('<h6>{}</h6>'.format(rs_line))
    for i, seq in enumerate(SeqIO.parse(open(f'{OUTPUT_DIR_T}/others/repeat_seq.fasta'), 'fasta')):
        if i == turnover:
            f.write('<h6>{}</h6>'.format(rs_line))
            f.write('<h6>{}</h6>'.format(rs_allele2_label))
            f.write('<h6>{}</h6>'.format(rs_line))
        if len(seq.seq) < 25000 and i+1 in subsampled_readnum:
            sumrep = summarize_tandem_repeats(str(seq.seq).replace('-', ''))
            f.write('<h6>{}</h6>'.format(seq.id + '\t' + colorize_repeats(sumrep)))
    f.write('<h1>4. Methylation profile</h1>')
    f.write('<h6>Repeat region = {}</h6>'.format(target_region))
    f.write(methylation_plot)
    f.write('<h4>Methylation key regions profile</h4>')
    f.write('<h6>{} in order.</h6>'.format(labels))
    f.write(methylation_key_plot)
    if meth_key_df is not None:
        f.write('<h5>Methylation key regions profile table</h5>')
        f.write(meth_key_df.to_html())
    f.write('</div>')
    f.write('</body>')
    f.write('</html>')

print(f'Data report is generated successfully. Please check the result.html file in {OUTPUT_DIR_T} directory.')