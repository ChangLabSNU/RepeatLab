configfile: 'config.yml'

SAMPLES = list(config['data']['sources'])
REFERENCE = config['options']['reference']

rule all:
    input:
        expand('basecalls/guppy/sup/{sample}/pass.fastq.gz', sample=SAMPLES),
        expand('alignments/{sample}_sup.sorted.bam', sample=SAMPLES),
        expand('alignments/{sample}_sup.sorted.bam.stats', sample=SAMPLES)
        #expand('logbam/RepBAM_{gene}.gmm_GapCorrection1_FlankLength30_SplitAndReAlign1_2_7_4_80_10_100_hg38_comp_{sample}_I0.120_D0.020_S0.020.log', sample=SAMPLES, gene=TARGETGENE),
        #expand('{sample}_{gene}_RepeatCount.png', sample=SAMPLES, gene=TARGETGENE)
        
rule guppy_basecall:
    output: 'basecalls/guppy/sup/{sample}/pass.fastq.gz'
    params:
        outputdir='basecalls/guppy/sup/{sample}/'
    threads: 10
    run:
       srcdir = config['data']['sources'][wildcards.sample]
       bcopts = config['options']['basecalling']
       shell(f'{config["programs"]["guppy"]} --recursive \
                   -x {bcopts["cuda_devices"]} \
                   -i {srcdir} -c {bcopts["model"]} \
                   -s {params.outputdir} --compress_fastq --disable_pings \
                   --gpu_runners_per_device {bcopts["gpu_runners"]} \
                   --num_callers {bcopts["num_callers"]}')
       shell('zcat {params.outputdir}/pass/fastq_*.fastq.gz | \
              bgzip -@ {threads} -c > {output}')

rule align:
    input: 
        fastq='basecalls/guppy/sup/{sample}/pass.fastq.gz',
        index=f'{REFERENCE}.genome.mm2.idx'
    output: temp('alignments/{sample}_sup.unsorted.bam')
    threads: 10
    shell: 'minimap2 -a -x map-ont -t {threads} --MD {input.index} {input.fastq} | \
            samtools view -b > {output}'
    
rule sort_bam:
    input: 'alignments/{sample}_sup.unsorted.bam'
    output:
        'alignments/{sample}_sup.sorted.bam',
        'alignments/{sample}_sup.sorted.bam.csi'
    threads: 10
    shell: 'samtools sort --write-index -@ {threads} -o {output} {input}'
    
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
        sampleid = SAMPLES
        shell('./cdrun.sh repeathmmenv {REPEATHMM_CMD} BAMinput --repeatName {targetgene} \
                --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID {sampleid} \
                --Onebamfile {input.bam} --outFolder {output} \
                --Patternfile {input.patternfile} \
                --hgfile {input.reference}')
                
rule plot_repeat_count_histogram:
    input: 'logbam/RepBAM_{gene}.gmm_GapCorrection1_FlankLength30_SplitAndReAlign1_2_7_4_80_10_100_hg38_comp_{sample}_I0.120_D0.020_S0.020.log'
    output: '{sample}_{gene}_RepeatCount.png'
    shell: 'python pyscripts/plotRepeatCount.py {input}'
    