
#################################### Defined by users #################################
configfile: "config/config_paired.yaml"    # Sets path to the config file
#######################################################################################

shell.prefix('set -euo pipefail; ')
shell.executable('/bin/bash')


THREADS=config['THREADS']

rule all: 
    input: 
        expand("star_output/featurecounts_{sample}.tsv", sample=config['FASTQ_PREFIX']),
        expand("hisat2_output/featurecounts_{sample}.tsv", sample=config['FASTQ_PREFIX'])


rule get_fastq:   # Creates fastq.gz files in fastq directory
    """
    This rule downloads SRA and converts to FASTQ files
    """
    output:
        expand("fastq/{out}.fastq.gz", out=config['FASTQ_LIST'])
    params:
        dic=config['SAMPLE'],
        reads=config['END'],
        sra=config['SRA']
    run:
        shell("fastq-dump --split-files {params.sra} --gzip -X 100000")    # with or without -X 
        for key, value in params.dic.items(): 
              for read in params.reads: 
                  shell("mv {key}_{read}.fastq.gz fastq/{value}_{read}.fastq.gz") 


rule align_star:   # Creates bam files in star_output directory"
    """
    This rule aligns the reads using STAR two-pass mode
    """
    input:
        config["GTF"],
        expand("fastq/{out}.fastq.gz", out=config['FASTQ_LIST'])
    output:
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", sample=config['FASTQ_PREFIX'])  
    params:
        indexing=config["INDEX_STAR"],
        indir=config['FASTQ_DIR'],
        files=config["FASTQ_PREFIX"], 
        read_ends=config['END'],
        ext=config['FASTQ_EXT']
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= params.indir + params.files[i] + "_1" + params.ext + " " 
            r2=""
            if len(params.read_ends) == 2: 
                r2= params.indir + params.files[i] + "_2" + params.ext + " " 
            shell("STAR --runThreadN {THREADS} "  
                    "--runMode alignReads "  
                    "--readFilesCommand zcat "
                    "--genomeDir {params.indexing} " 
                    "--sjdbGTFfile {input[0]} "  
                    "--sjdbOverhang 100 "  
                    "--readFilesIn {r1}{r2}"  
                    "--outFileNamePrefix star_output/{p} "
                    "--outFilterType BySJout "  
                    "--outFilterMultimapNmax 20 "
                    "--alignSJoverhangMin 8 "
                    "--alignSJDBoverhangMin 1 "
                    "--outFilterMismatchNmax 999 "
                    "--outFilterMismatchNoverReadLmax 0.04 "
                    "--alignIntronMin 20 "
                    "--alignIntronMax 1000000 "
                    "--outSAMunmapped None "
                    "--outSAMtype BAM "
                    "SortedByCoordinate "
                    "--quantMode GeneCounts "
                    "--twopassMode Basic "
                    "--chimOutType Junctions")     


rule align_hisat2:    # Creates bam files in hisat2_output directory"
    """
    This rule aligns the reads using HISAT2    
    """
    input:
        expand("fastq/{out}.fastq.gz", out=config['FASTQ_LIST'])
    output:
        expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX']),
        temp(expand("hisat2_output/{sample}.sam", sample=config['FASTQ_PREFIX']))
    params:
        indir=config['FASTQ_DIR'],
        files=config["FASTQ_PREFIX"], 
        read_ends=config['END'],
        ext=config['FASTQ_EXT'],
        indexing=config['INDEX_HISAT']
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= params.indir + params.files[i] + "_1" + params.ext 
            r2=""
            read="-U " + r1
            if len(params.read_ends) == 2: 
                r2= params.indir + params.files[i] + "_2" + params.ext  
                read="-1 " + r1 + " -2 " + r2
            shell("hisat2 -q -p {THREADS} "
                  "--seed 23 "
                  "--summary-file hisat2_output/summary_{p}.txt "
                  "-x {params.indexing} "
                  "{read} "
                  "-S hisat2_output/{p}.sam && "
                  "samtools view -bS "
                  "-@ {THREADS} "
                  "hisat2_output/{p}.sam > hisat2_output/{p}.bam")

rule featurecounts:
    """
    This rule assesses read counts using featureCounts
    """
    input:
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", sample=config['FASTQ_PREFIX']), 
        expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX'])
    params:
        ends=config['END'],
        gtf=config['GTF']
    output:
        expand("star_output/featurecounts_{sample}.tsv", sample=config['FASTQ_PREFIX']),
        expand("hisat2_output/featurecounts_{sample}.tsv", sample=config['FASTQ_PREFIX'])
    run:
        end=""
        if len(params.ends) == 2:
            end='-p'
        for i in range(len(output)):
            bam=input[i]
            counts=output[i]
            # Additional featureCounts flags
            # -F: format of the annotation file. 'GTF' by default.
            # -g: attribute type. 'gene_id' by default. 
            # -L: set for long-read inputs 
            # -s: strand-specificity. 0 (unstranded & default), 1 (stranded), 2 (reversely stranded)
            shell("featureCounts {end} "    
                  "-s0 "
                  "-T {THREADS} "
                  "-a {params.gtf} "
                  "-o {counts} "
                  "{bam} &> {counts}.log")


    
    
