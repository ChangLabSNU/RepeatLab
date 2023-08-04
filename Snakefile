configfile: 'config.yml'

SAMPLES = list(config['data']['sources'])
SAMPLES_MTPLX = list(config['data']['sources_multiplex'])
REFERENCE = config['options']['reference']

rule all:
    input:
        expand('raw_pod5/{sample}/', sample=SAMPLES+SAMPLES_MTPLX),
        expand('basecalls/dorado/fast/{sample}/{sample}.fast-called.fastq.gz', sample=SAMPLES+SAMPLES_MTPLX),
        expand('basecalls/dorado/fast/{sample_mtplx}/{sample_mtplx}.demultiplexed.fast-called.fastq.gz', sample_mtplx=SAMPLES_MTPLX),
        expand('alignments/{sample}.fast-aligned.sorted.bam', sample=SAMPLES+SAMPLES_MTPLX),
        expand('alignments/{sample}.fast-aligned.sorted.bam.bai', sample=SAMPLES+SAMPLES_MTPLX),
        expand('alignments/{sample}.fast-aligned.sorted.bam.stats', sample=SAMPLES+SAMPLES_MTPLX),
        expand('on-target/readID_list/{analysis}.ontarget.readID.txt', analysis=config['analysis']),
        expand('basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam', analysis=config['analysis']),
        expand('basecalls/dorado/sup_v3.6/{analysis}.sup-called.summary.txt', analysis=config['analysis']),
        expand('alignments/{analysis}.sup-aligned.sorted.bam', analysis=config['analysis']),
        expand('alignments/{analysis}.sup-aligned.sorted.bam.bai', analysis=config['analysis']),
        expand('analyses/{analysis}/result.log', analysis=config['analysis']),
        expand('on-target/phased_readID_list/{analysis}.allele{num}.readID.txt', analysis=config['analysis'], num=[1,2])

rule convert_fast5_to_pod5:
    output: directory('raw_pod5/{sample}/')
    priority: 100
    run:
        if wildcards.sample in SAMPLES:
            srcdir = config['data']['sources'][wildcards.sample]
        elif wildcards.sample in SAMPLES_MTPLX:
            srcdir = config['data']['sources_multiplex'][wildcards.sample]
        shell(f'conda run --no-capture-output -n {config["programs"]["pod5_condaenv"]} \
                pod5 convert fast5 -o {output} -r {srcdir} --one-to-one {srcdir}')

ruleorder: dorado_basecall_first > demultiplex

rule dorado_basecall_first:
    input: 'raw_pod5/{sample}/'
    output: 'basecalls/dorado/fast/{sample}/{sample}.fast-called.fastq.gz'
    threads: 10
    priority: 99
    run:
        bcopts = config['options']['first_basecalling']
        shell(f'conda run --no-capture-output -n {config["programs"]["dorado_condaenv"]} \
                dorado basecaller -x {bcopts["cuda_devices"]} -r --emit-fastq {bcopts["model"]} {input} \
                | bgzip -@ {threads} -c > {output}')

rule demultiplex:
    output: 'basecalls/dorado/fast/{sample_mtplx}/{sample_mtplx}.demultiplexed.fast-called.fastq.gz'
    threads: 10
    priority: 98
    run:
        inputdir='basecalls/dorado/fast/{wildcards.sample_mtplx}/'
        outputdir='basecalls/dorado/fast/{wildcards.sample_mtplx}_multiplexing/'
        bcopts = config['options']['demultiplexing']
        barcode_kit = bcopts['barcode_kit']
        barcode_num = bcopts['barcode_num'][wildcards.sample_mtplx]
        shell(f'{config["programs"]["guppy_barcoder"]} -i {inputdir} --disable_pings -t {threads} -x cuda:6,7 \
                --barcode_kits {barcode_kit} --compress_fastq -s {outputdir}')
        shell(f'zcat {outputdir}/{barcode_num}/fastq_*.fastq.gz | \
                bgzip -@ {threads} -c > {output}')     

rule align_first:
    input:
        index=f'{REFERENCE}.genome.mm2.idx'
    output: temp('alignments/{sample}.fast-aligned.unsorted.bam')
    threads: 32
    priority: 97
    run:
        if wildcards.sample in SAMPLES:
            fastq = 'basecalls/dorado/fast/{wildcards.sample}/{wildcards.sample}.fast-called.fastq.gz'
        elif wildcards.sample in SAMPLES_MTPLX:
            fastq = 'basecalls/dorado/fast/{wildcards.sample}/{wildcards.sample}.demultiplexed.fast-called.fastq.gz'
        shell(f'conda run --no-capture-output -n {config["programs"]["dorado_condaenv"]} \
                dorado aligner -t {threads} {input.index} {fastq} > {output}')

rule sort_alignment_first:
    input: 'alignments/{sample}.fast-aligned.unsorted.bam'
    output: 'alignments/{sample}.fast-aligned.sorted.bam'
    threads: 10
    priority: 96
    shell: 'samtools sort -@ {threads} -o {output} {input}'

rule index_alignment_first:
    input: 'alignments/{sample}.fast-aligned.sorted.bam'
    output: 'alignments/{sample}.fast-aligned.sorted.bam.bai'
    priority: 95
    shell: 'samtools index {input}'

rule stats_alignment_first:
    input: 'alignments/{sample}.fast-aligned.sorted.bam'
    output: 'alignments/{sample}.fast-aligned.sorted.bam.stats'
    priority: 94
    shell: 'samtools stats {input} > {output}'

def prepare_inputbam_ontarget(wildcards):
    analysis_setting = config['analysis'][wildcards.analysis]
    sample = analysis_setting['source']
    return {
        'pod5': f'raw_pod5/{sample}/',
        'bam': f'alignments/{sample}.fast-aligned.sorted.bam'
    }

rule extract_ontarget_reads:
    input: unpack(prepare_inputbam_ontarget)
    output: 'on-target/readID_list/{analysis}.ontarget.readID.txt'
    priority: 93
    run:
        target_gene = str(wildcards.analysis).split('-')[1]
        target_region = config['options']['on-target_extraction'][target_gene]['region']
        shell('samtools view {input.bam} {target_region} | cut -f1 | sort | uniq > {output}')

rule dorado_basecall_second:
    input: unpack(prepare_inputbam_ontarget)
    output: 'basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam'
    threads: 10
    priority: 92
    run:
        readID = 'on-target/readID_list/{wildcards.analysis}.ontarget.readID.txt'
        bcopts = config['options']['second_basecalling']
        shell(f'conda run --no-capture-output -n {config["programs"]["dorado_condaenv"]} \
                dorado basecaller -x {bcopts["cuda_devices"]} --batchsize 128 --chunksize 80000 \
                 -r --emit-sam --emit-moves {bcopts["model"]} {input.pod5} -l {readID} \
                | samtools view -b -o {output}')

rule dorado_summary:
    input: 'basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam'
    output: 'basecalls/dorado/sup_v3.6/{analysis}.sup-called.summary.txt'
    priority: 91
    run:
        shell(f'conda run --no-capture-output -n {config["programs"]["dorado_condaenv"]} \
                dorado summary {input} > {output}')

rule align_second:
    input:
        bam='basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam',
        index=f'{REFERENCE}.genome.mm2.idx'
    output: temp('alignments/{analysis}.sup-aligned.unsorted.bam')
    threads: 32
    priority: 90
    run:
        shell(f'conda run --no-capture-output -n {config["programs"]["dorado_condaenv"]} \
                dorado aligner -t {threads} {input.index} {input.bam} > {output}')

rule sort_alignment_second:
    input: 'alignments/{analysis}.sup-aligned.unsorted.bam'
    output: 'alignments/{analysis}.sup-aligned.sorted.bam'
    threads: 10
    priority: 89
    shell: 'samtools sort -@ {threads} -o {output} {input}'

rule index_alignment_second:
    input: 'alignments/{analysis}.sup-aligned.sorted.bam'
    output: 'alignments/{analysis}.sup-aligned.sorted.bam.bai'
    priority: 88
    shell: 'samtools index {input}'

def prepare_inputs_repeathmm(wildcards):
    analysis_setting = config['analysis'][wildcards.analysis]
    sample = analysis_setting['source']
    return {
        'patternfile': config['options']['repeat_presets'],
        'bam': f'alignments/{sample}_sup.sorted.bam',
        'bai': f'alignments/{sample}_sup.sorted.bam.bai',
        'reference': os.path.abspath(f'{REFERENCE}.genome.fa.gz')
    }

rule run_repeathmm:
    input:
        patternfile=config['options']['repeat_presets'],
        bam='alignments/{analysis}.sup-aligned.sorted.bam',
        bai='alignments/{analysis}.sup-aligned.sorted.bam.bai',
        reference=os.path.abspath(f'{REFERENCE}.genome.fa.gz')
    output: 'analyses/{analysis}/result.log'
    params:
        tmpdir='analyses/{analysis}/tmp',
        logbam='logbam'
    priority: 87
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
                python2.7 {config["programs"]["repeathmm"]} BAMinput --repeatName {gene} \
                --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID {wildcards.analysis} \
                --Onebamfile {input.bam} --outFolder {params.tmpdir}/ \
                --Patternfile {input.patternfile} \
                --hgfile {input.reference} \
                --SplitAndReAlign 2 --SeqTech Nanopore --outlog DEBUG')
        shell('mv {params.logbam}/RepBAM_{gene}.gmm*_{sample}_*.log {output}')

rule phasing:
    output: 'on-target/phased_readID_list/{analysis}.allele{num}.readID.txt'
    run: 
        shell(f'python scripts/phasing.py {wildcards.analysis}')

