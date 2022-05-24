configfile: 'config.yml'

SAMPLES = list(config['data']['sources'])
REFERENCE = config['options']['reference']

rule all:
    input:
        expand('basecalls/guppy/sup/{sample}/pass.fastq.gz', sample=SAMPLES),
        expand('alignments/{sample}_sup.sorted.bam', sample=SAMPLES),
        expand('alignments/{sample}_sup.sorted.bam.stats', sample=SAMPLES),
        expand('analyses/{analysis}/result.log', analysis=config['analysis']),
        expand('plots/RepeatCount-{analysis}.png', analysis=config['analysis'])

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
                   --min_qscore {bcopts["min_qscore"]} \
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

def prepare_inputs_repeathmm(wildcards):
    analysis_setting = config['analysis'][wildcards.analysis]
    sample = analysis_setting['source']
    return {
        'patternfile': config['options']['repeat_presets'],
        'bam': f'alignments/{sample}_sup.sorted.bam',
        'bai': f'alignments/{sample}_sup.sorted.bam.bai',
        'reference': os.path.abspath(f'{REFERENCE}.genome.fa.gz')
    }

rule count_repeat_repeathmm:
    input: unpack(prepare_inputs_repeathmm)
    output: 'analyses/{analysis}/result.log'
    params:
        tmpdir='analyses/{analysis}/tmp',
        logbam='logbam'
    run:
        analysis_setting = config['analysis'][wildcards.analysis]
        sample = analysis_setting['source']
        gene = analysis_setting['target']
        if not os.path.isdir(params.tmpdir):
            os.makedirs(params.tmpdir)
        if os.path.isdir(params.logbam):
            shell('rm -rf {params.logbam}')
        os.makedirs(params.logbam)

        shell(f'conda run --no-capture-output -n {config["programs"]["repeathmm_condaenv"]} \
                python2 {config["programs"]["repeathmm"]} BAMinput --repeatName {gene} \
                --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID {wildcards.analysis} \
                --Onebamfile {input.bam} --outFolder {params.tmpdir}/ \
                --Patternfile {input.patternfile} \
                --hgfile {input.reference} \
                --SplitAndReAlign 2 --SeqTech Nanopore --outlog DEBUG')
        shell('mv {params.logbam}/RepBAM_{gene}.gmm*_{sample}_*.log {output}')

rule plot_repeat_count_histogram:
    input: 'analyses/{analysis}/result.log'
    output: 'plots/RepeatCount-{analysis}.png'
    shell: 'python pyscripts/plotRepeatCount.py {input} {output}'
