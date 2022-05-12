SOURCES = {
    'DM12' : '20220203-SSLabMed-Nanopore/fast5/F_DMPK_DM12/'
}

SAMPLE = ['DM12']

TARGETGENE = ['DMPK']

REFERENCE = 'GRCh38.v39'

CONFIG = 'dna_r9.4.1_450bps_sup.cfg'

GUPPY_CMD = '~/gh/ont-guppy/bin/guppy_basecaller'
REPEATHMM_CMD = 'python ~/gh/RepeatHMM/bin/repeatHMM.py'

rule all:
    input:
        expand('basecall/guppy/sup/{sample}/basecalls.fastq.gz', sample=SAMPLE),
        expand('alignments/{sample}_sup.unsorted.bam', sample=SAMPLE),
        expand('alignments/{sample}_sup.sorted.bam', sample=SAMPLE),
        expand('alignments/{sample}_sup.sorted.bam.bai', sample=SAMPLE),
        expand('alignments/{sample}_sup.sorted.bam.stats', sample=SAMPLE),
        expand('logbam/RepBAM_{gene}.gmm_GapCorrection1_FlankLength30_SplitAndReAlign1_2_7_4_80_10_100_hg38_comp_{sample}_I0.120_D0.020_S0.020.log', sample=SAMPLE, gene=TARGETGENE),
        expand('{sample}_{gene}_RepeatCount.png', sample=SAMPLE, gene=TARGETGENE)
        
rule guppy_basecall:
    output: 'basecall/guppy/sup/{sample}/basecalls.fastq.gz'
    params:
        outputdir='basecall/guppy/sup/{sample}/'
    run:
       srcdir = SOURCES[wildcards.sample]
       config = CONFIG
       gpus = 'cuda:6,7'
       shell('{GUPPY_CMD} --recursive \
                   -x {gpus} -i {srcdir} -c {config} \
                   -s {params.outputdir} --compress_fastq --disable_pings \
                   --gpu_runners_per_device 2 --num_callers 2')
       shell('zcat {params.outputdir}/pass/fastq_*.fastq.gz | bgzip -@ 10 -c /dev/stdin > {output}')

rule align:
    input: 
        fastq = 'basecall/guppy/sup/{sample}/basecalls.fastq.gz',
        index = 'reference/GRCh38.v39.genome.mm2.idx'
    output: 'alignments/{sample}_sup.unsorted.bam'
    threads: 10
    shell: 'minimap2 -a -x map-ont -t {threads} --MD {input.index} {input.fastq} | \
            samtools view -b > {output}'
    
rule sort_bam:
    input: 'alignments/{sample}_sup.unsorted.bam'
    output: 'alignments/{sample}_sup.sorted.bam'
    threads: 10
    shell: 'samtools sort -@ {threads} -o {output} {input}'
    
rule index_bam:
    input: 'alignments/{sample}_sup.sorted.bam'
    output: 'alignments/{sample}_sup.sorted.bam.bai'
    shell: 'samtools index {input}'
    
rule stats_bam:
    input: 'alignments/{sample}_sup.sorted.bam'
    output: 'alignments/{sample}_sup.sorted.bam.stats'
    shell: 'samtools stats {input} > {output}'
    
rule count_repeat_repeathmm:
    input:
        patternfile = 'GRCh38.v39.predefined.pa',
        bam = 'alignments/{sample}_sup.sorted.bam',
        reference = 'reference/GRCh38.v39.genome.fa.gz'
    output: 'logbam/RepBAM_{gene}.gmm_GapCorrection1_FlankLength30_SplitAndReAlign1_2_7_4_80_10_100_hg38_comp_{sample}_I0.120_D0.020_S0.020.log'
    run: 
        targetgene = TARGETGENE
        sampleid = SAMPLE
        shell('./cdrun.sh repeathmmenv {REPEATHMM_CMD} BAMinput --repeatName {targetgene} \
                --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID {sampleid} \
                --Onebamfile {input.bam} --outFolder {output} \
                --Patternfile {input.patternfile} \
                --hgfile {input.reference}')
                
rule plot_repeat_count_histogram:
    input: 'logbam/RepBAM_{gene}.gmm_GapCorrection1_FlankLength30_SplitAndReAlign1_2_7_4_80_10_100_hg38_comp_{sample}_I0.120_D0.020_S0.020.log'
    output: '{sample}_{gene}_RepeatCount.png'
    shell: 'python plotRepeatCount.py {input}'
    