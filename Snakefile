import os
import glob

configfile: 'config.yml'

SAMPLES = list(config['data'].get('sources', {}) or {})
SAMPLES_MTPLX = list(config['data'].get('sources_multiplex', {}) or {})
ALL_SAMPLES = SAMPLES + SAMPLES_MTPLX
REFERENCE = config['options']['reference']

if len(ALL_SAMPLES) == 0:
    raise ValueError("ERROR: No valid samples found in 'sources' or 'sources_multiplex' in config.yml.")

rule all:
    input:
        expand('raw_pod5/{sample}/', sample=ALL_SAMPLES),
        expand('basecalls/dorado/fast/{sample}/{sample}.fast-called.bam', sample=ALL_SAMPLES),
        (expand('basecalls/dorado/fast/{sample_mtplx}/{sample_mtplx}.fast-called.bam', sample_mtplx=SAMPLES_MTPLX) if SAMPLES_MTPLX else []),
        expand('basecalls/dorado/fast/{sample}/{sample}.fast-called.summary.txt', sample=ALL_SAMPLES),
        expand('alignments/{sample}.fast-aligned.sorted.bam', sample=ALL_SAMPLES),
        expand('alignments/{sample}.fast-aligned.sorted.bam.bai', sample=ALL_SAMPLES),
        expand('alignments/{sample}.fast-aligned.sorted.bam.stats', sample=ALL_SAMPLES),
        expand('on-target/readID_list/{analysis}.ontarget.readID.txt', analysis=config['analysis']),
        expand('basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam', analysis=config['analysis']),
        expand('basecalls/dorado/sup_v3.6/{analysis}.sup-called.summary.txt', analysis=config['analysis']),
        expand('alignments/{analysis}.sup-aligned.sorted.bam', analysis=config['analysis']),
        expand('alignments/{analysis}.sup-aligned.sorted.bam.bai', analysis=config['analysis']),
        expand('analyses/{analysis}/result.log', analysis=config['analysis']),
        expand('on-target/phased_readID_list/{analysis}.allele{num}.readID.txt', analysis=config['analysis'], num=[1,2]),
        expand('reference/{analysis}_target_ref_genome.{suffix}', suffix=['fa', 'mmi'], analysis=config['analysis']),
        expand('on-target/phased_readID_list/{analysis}.allele{num}.readID.txt', analysis=config['analysis'], num=[1,2]),
        expand('basecalls/mod/{analysis}.allele{num}.meth-called.bam', analysis=config['analysis'], num=[1,2]),
        expand('alignments/{analysis}.allele{num}.meth-aligned.bam', analysis=config['analysis'], num=[1,2]),
        expand('alignments/{analysis}.allele{num}.meth-aligned.sorted.bam', analysis=config['analysis'], num=[1,2]),
        expand('meth-profiles/{analysis}.allele{num}.meth-called.bedMethyl', analysis=config['analysis'], num=[1,2]),
        expand('meth-profiles/{analysis}.allele{num}.meth-called.tsv', analysis=config['analysis'], num=[1,2]),        
        expand('analyses/{analysis}/result.html', analysis=config['analysis'])

rule convert_fast5_to_pod5:
    input: 
        lambda wildcards: glob.glob(f"{config['data']['sources'][wildcards.sample]}/*.fast5")
    output: directory('raw_pod5/{sample}/')
    priority: 100
    run:
        srcdir = config['data']['sources'][wildcards.sample]
        fast5_files = glob.glob(f'{srcdir}/*.fast5')
        if not fast5_files:
            # Skip conversion if no .fast5 files are found
            shell(f'mkdir -p {output}')
            pod5_files = glob.glob(f'{srcdir}/*.pod5')
            if pod5_files:
                shell(f'cp -r {srcdir}/*.pod5 {output}')
            else:
                raise FileNotFoundError(f'No .fast5 or .pod5 files found in {srcdir}')
        else:
            shell(f'conda run --no-capture-output -n {config["programs"]["pod5_condaenv"]} \
                    pod5 convert fast5 -o {output} -r {srcdir} --one-to-one {srcdir}')

if SAMPLES_MTPLX:
    ruleorder: dorado_basecall_first > demultiplex

rule dorado_basecall_first:
    input: 'raw_pod5/{sample}/'
    output: 'basecalls/dorado/fast/{sample}/{sample}.fast-called.bam'
    threads: 10
    priority: 99
    run:
        bcopts = config['options']['first_basecalling']
        shell(f'{config["programs"]["dorado"]} basecaller -x {bcopts["cuda_devices"]} -r --emit-sam {bcopts["model"]} {input} \
                | samtools view -b -o {output}')

if SAMPLES_MTPLX:
    rule demultiplex:
        input: 'basecalls/dorado/fast/{sample_mtplx}/{sample_mtplx}.fast-called.bam'
        output: 'basecalls/dorado/fast/{sample_mtplx}/{sample_mtplx}.fast-called.bam'
        threads: 10
        priority: 98
        run:
            demuxdir='basecalls/dorado/fast/demux/'
            bcopts = config['options']['demultiplexing']
            barcode_kit = bcopts['barcode_kit']
            barcode_num = bcopts['barcode_num'][wildcards.sample_mtplx]

            if not os.path.exists(demuxdir):
                shell(f'conda run --no-capture-output -n {config["programs"]["dorado_condaenv"]} \
                        dorado demux --recursive --kit-name {barcode_kit} --output-dir {demuxdir} {input}')
            shell(f'cp {demuxdir}/{barcode_kit}_{barcode_num}.bam {output}') 

rule dorado_summary_first:
    input: 'basecalls/dorado/fast/{sample}/{sample}.fast-called.bam'
    output: 'basecalls/dorado/fast/{sample}/{sample}.fast-called.summary.txt'
    priority: 97
    run:
        shell(f'{config["programs"]["dorado"]} summary {input} > {output}')

rule align_first:
    input:
        index=f'{REFERENCE}.genome.mm2.idx'
    output: temp('alignments/{sample}.fast-aligned.unsorted.bam')
    threads: 32
    priority: 96
    run:
        fast_basecalled = 'basecalls/dorado/fast/{wildcards.sample}/{wildcards.sample}.fast-called.bam'
        shell(f'{config["programs"]["dorado"]} aligner -t {threads} {input.index} {fast_basecalled} > {output}')

rule sort_alignment_first:
    input: 'alignments/{sample}.fast-aligned.unsorted.bam'
    output: 'alignments/{sample}.fast-aligned.sorted.bam'
    threads: 10
    priority: 95
    shell: 'samtools sort -@ {threads} -o {output} {input}'

rule index_alignment_first:
    input: 'alignments/{sample}.fast-aligned.sorted.bam'
    output: 'alignments/{sample}.fast-aligned.sorted.bam.bai'
    priority: 94
    shell: 'samtools index {input}'

rule stats_alignment_first:
    input: 'alignments/{sample}.fast-aligned.sorted.bam'
    output: 'alignments/{sample}.fast-aligned.sorted.bam.stats'
    priority: 93
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
    priority: 92
    run:
        target_gene = str(wildcards.analysis).split('-')[1]
        target_region = config['options']['on-target_extraction'][target_gene]['region']
        shell('samtools view {input.bam} {target_region} | cut -f1 | sort | uniq > {output}')

rule dorado_basecall_second:
    input: unpack(prepare_inputbam_ontarget)
    output: 'basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam'
    threads: 10
    priority: 91
    run:
        readID = 'on-target/readID_list/{wildcards.analysis}.ontarget.readID.txt'
        bcopts = config['options']['second_basecalling']
        shell(f'{config["programs"]["dorado"]} basecaller -x {bcopts["cuda_devices"]} --batchsize 128 --chunksize 60000 \
                 -r --emit-sam {bcopts["model"]} {input.pod5} -l {readID} \
                | samtools view -b -o {output}')

rule dorado_summary_second:
    input: 'basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam'
    output: 'basecalls/dorado/sup_v3.6/{analysis}.sup-called.summary.txt'
    priority: 90
    run:
        shell(f'{config["programs"]["dorado"]} summary {input} > {output}')

rule align_second:
    input:
        bam='basecalls/dorado/sup_v3.6/{analysis}.sup-called.bam',
        index=f'{REFERENCE}.genome.mm2.idx'
    output: temp('alignments/{analysis}.sup-aligned.unsorted.bam')
    threads: 32
    priority: 89
    run:
        shell(f'{config["programs"]["dorado"]} aligner -t {threads} {input.index} {input.bam} > {output}')

rule sort_alignment_second:
    input: 'alignments/{analysis}.sup-aligned.unsorted.bam'
    output: 'alignments/{analysis}.sup-aligned.sorted.bam'
    threads: 10
    priority: 88
    shell: 'samtools sort -@ {threads} -o {output} {input}'

rule index_alignment_second:
    input: 'alignments/{analysis}.sup-aligned.sorted.bam'
    output: 'alignments/{analysis}.sup-aligned.sorted.bam.bai'
    priority: 87
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
    priority: 86
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
                --MinSup 0 --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID {wildcards.analysis} \
                --Onebamfile {input.bam} --outFolder {params.tmpdir}/ \
                --Patternfile {input.patternfile} \
                --hgfile {input.reference} \
                --SplitAndReAlign 2 --SeqTech Nanopore --outlog DEBUG')
        shell('mv {params.logbam}/RepBAM_{gene}.gmm*_{sample}_*.log {output}')

rule phasing:
    output: 'on-target/phased_readID_list/{analysis}.allele{num}.readID.txt'
    priority: 85
    run: 
        shell(f'python scripts/phasing.py {wildcards.analysis}')

rule uncompress_genome:
    input: f'{REFERENCE}.genome.fa.gz'
    output: temp('reference/full_ref_genome.fa.gz')
    priority: 84
    shell: 'gunzip -c {input} > {output}'

rule extract_genome_chromosomes:
    input: 'reference/full_ref_genome.fa.gz',
    output: 'reference/{analysis}_target_ref_genome.fa'
    priority: 83
    run:
        analysis_setting = config['analysis'][wildcards.analysis]
        gene = analysis_setting['target']
        target_chromosome=config['options']['on-target_extraction'][gene]['chromosome']

        shell('samtools faidx {input} {target_chromosome} > {output}')

rule build_minimap_index:
    input: 'reference/{analysis}_target_ref_genome.fa'
    output: 'reference/{analysis}_target_ref_genome.mmi'
    priority: 82
    shell: 'minimap2 -x map-ont -k 13 -w 20 -d {output} {input}'

rule basecall_dorado_5mC:
    input:
        read_ids='on-target/phased_readID_list/{analysis}.allele{num}.readID.txt'
    output: 'basecalls/mod/{analysis}.allele{num}.meth-called.bam'
    priority: 81
    run:
        sample = str(wildcards.analysis).split('-')[0]
        datadir=f'raw_pod5/{sample}/'
        bcopts = config['options']['methylation_calling']
        shell(f'{config["programs"]["dorado"]} basecaller -x {bcopts["cuda_devices"]} --batchsize 128 --chunksize 10000 \
                -l {input.read_ids} -r \
                --modified-bases-models {bcopts["dorado_5mCG_model"]} \
                --emit-sam {bcopts["dorado_original_model"]} \
                {datadir} | samtools view -b -o {output}')

rule align_sequences:
    input: 
        basecalls='basecalls/mod/{analysis}.allele{num}.meth-called.bam',
        index='reference/{analysis}_target_ref_genome.mmi'
    output: 'alignments/{analysis}.allele{num}.meth-aligned.bam'
    priority: 80
    run: 
        threads=32
        shell(f'{config["programs"]["dorado"]} aligner -t {threads} {input.index} {input.basecalls} > {output}')

rule sort_alignments:
    input: 'alignments/{analysis}.allele{num}.meth-aligned.bam'
    output: 'alignments/{analysis}.allele{num}.meth-aligned.sorted.bam'
    priority: 79
    shell: 'samtools sort --write-index -o {output} {input}'

rule pileup_methylated_sites:
    input:
        alignments='alignments/{analysis}.allele{num}.meth-aligned.sorted.bam',
        reference='reference/{analysis}_target_ref_genome.fa'
    output: 'meth-profiles/{analysis}.allele{num}.meth-called.bedMethyl'
    priority: 78
    run: 
        shell(f'{config["programs"]["modkit"]} pileup --no-filtering --combine-strands --cpg \
                --ref {input.reference} {input.alignments} {output}')

rule tabularize_bedmethyl:
    input: 'meth-profiles/{analysis}.allele{num}.meth-called.bedMethyl'
    output: 'meth-profiles/{analysis}.allele{num}.meth-called.tsv'
    priority: 77
    shell: '(echo "chrom\tposition\tcoverage\tmethylation_pct"; \
             awk -vOFS=\'\\t\' \'{{print $1, $2, $10, $11}}\' \
             {input}) > {output}'

rule result_visualize:
    input: 'analyses/{analysis}/result.log'
    output: 'analyses/{analysis}/result.html'
    priority: 76
    run:
        shell('python scripts/visualization.py {wildcards.analysis}')
