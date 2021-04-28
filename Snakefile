
#################################### Defined by users #################################
configfile: "config/config_paired.yaml"    # Sets path to the config file
#######################################################################################

shell.prefix('set -euo pipefail; ')
shell.executable('/bin/bash')


THREADS=config['THREADS']

rule all: 
    input: 
        "star_output/featurecounts.tsv",
        "hisat2_output/featurecounts.tsv"

rule get_fastq:   # Creates fastq.gz files in fastq directory
    """
    This rule downloads SRA and converts to FASTQ files
    """
    output:
        expand("fastq/{out}.fastq.gz", out=config['FASTQ_LIST'])  # Gzipped FASTQ files from SRA 
    params:
        dic=config['SAMPLE'],   # Sample dictionary
        reads=config['END'],    # Reads (e.g. [1] or [1, 2]) 
        sra=config['SRA']       # SRA number 
    run:
        shell("fastq-dump --split-files {params.sra} --gzip -X 100000")    # with or without -X 
        for key, value in params.dic.items(): 
              for read in params.reads: 
                  shell("mv {key}_{read}.fastq.gz fastq/{value}_{read}.fastq.gz") 


rule get_reference:
    params:
        gen_link=config['REFERENCE_LINK']['GENOME'][0],   # Gencode reference genome file link 
        gen_name=config['REFERENCE_LINK']['GENOME'][1],   # Output reference genome location & name 
        anno_link=config['REFERENCE_LINK']['ANNOTATION'][0],  # Gencode GTF (annotation) file link
        anno_name=config['REFERENCE_LINK']['ANNOTATION'][1]   # Output GTF file location & name
    output:
        gen=expand("reference/{gen}", gen=config['REFERENCE_LINK']['GENOME'][2]),  # Decompressed reference genome file 
        anno=expand("reference/{anno}", anno=config['REFERENCE_LINK']['ANNOTATION'][2])  # Decompressed GTF file  
    shell:
        "set +o pipefail; "
        "wget -c {params.gen_link} -O reference/{params.gen_name} && "
        "wget -c {params.anno_link} -O reference/{params.anno_name} && "
        "gzip -d reference/*.gz"

rule index_hisat2:
    input: 
        expand("reference/{gen}", gen=config['REFERENCE_LINK']['GENOME'][2])    # Decompressed reference genome file
    output:
        expand("reference/hisat2_index/hisat2_index.{number}.ht2", number=[x+1 for x in range(8)])   # HISAT2 indexing files
    shell:
        "set +o pipefail; "
        "hisat2-build -f -o 4 "
        "-p {THREADS} "
        "--seed 67 "
        "{input} "
        "hisat2_index && "
        "mv *.ht2 reference/hisat2_index"


rule index_star:
    input:
        fa=expand("reference/{gen}", gen=config['REFERENCE_LINK']['GENOME'][2]),  # Decompressed reference genome file
        gtf=expand("reference/{anno}", anno=config['REFERENCE_LINK']['ANNOTATION'][2])  # Decompressed GTF file
    output:
        "reference/star_index/Genome",   # STAR indexing files
        "reference/star_index/SA",       # STAR indexing files
        "reference/star_index/SAindex"   # STAR indexing files
    shell:
        "set +o pipefail; "
        "STAR --runThreadN {THREADS} "
        "--runMode genomeGenerate "
        "--genomeDir reference/star_index "
        "--genomeFastaFiles {input.fa} "
        "--sjdbGTFfile {input.gtf}"




rule align_hisat2:    # Creates bam files in hisat2_output directory"
    """
    This rule aligns the reads using HISAT2    
    """
    input:
        expand("fastq/{out}.fastq.gz", out=config['FASTQ_LIST']),  # Gzipped FASTQ files
        expand("reference/hisat2_index/hisat2_index.{number}.ht2", number=[x+1 for x in range(8)])  # HISAT2 indexing files
    output:
        expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX']),   # Bam files
        temp(expand("hisat2_output/{sample}.sam", sample=config['FASTQ_PREFIX']))
    params:
        files=config["FASTQ_PREFIX"],  # e.g. Ctrl, Treatment
        read_ends=config['END'],       # e.g. [1] or [1, 2]
        ext=config['FASTQ_EXT'],       # extension of the FASTQ files (e.g. .fastq.gz)
        indexing=config['INDEX_HISAT'] # HISAT2 indexing location and file name prefix
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= "fastq/" + params.files[i] + "_1" + params.ext 
            r2=""
            read="-U " + r1
            if len(params.read_ends) == 2: 
                r2= "fastq/" + params.files[i] + "_2" + params.ext  
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



rule align_star:   # Creates bam files in star_output directory"
    """
    This rule aligns the reads using STAR two-pass mode
    """
    input:
        expand("reference/{gen}", gen=config['REFERENCE_LINK']['ANNOTATION'][2]),  # Decompressed GTF file
        expand("fastq/{out}.fastq.gz", out=config['FASTQ_LIST']),                  # Gzipped FASTQ files
        "reference/star_index/Genome",
        "reference/star_index/SA",
        "reference/star_index/SAindex"
    output:
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", sample=config['FASTQ_PREFIX'])  # Bam files
    params:
        indexing=config["INDEX_STAR"],  # STAR indexing file directory
        files=config["FASTQ_PREFIX"],   # e.g. Ctrl, Treatment
        read_ends=config['END'],        # e.g. [1] or [1, 2]
        ext=config['FASTQ_EXT']         # extension of the FASTQ files (e.g. fastq.gz)
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= "fastq/" + params.files[i] + "_1" + params.ext + " " 
            r2=""
            if len(params.read_ends) == 2: 
                r2= "fastq/" + params.files[i] + "_2" + params.ext + " " 
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



rule featurecounts:
    """
    This rule assesses read counts using featureCounts
    """
    input:
        star=expand("star_output/{sample}Aligned.sortedByCoord.out.bam", sample=config['FASTQ_PREFIX']),  # Bam files from STAR
        hisat2=expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX'])    # Bam files from HISAT2
    params:
        ends=config['END'], # Read ends (e.g. [1] or [1, 2])
        gtf=expand("reference/{anno}", anno=config['REFERENCE_LINK']['ANNOTATION'][2])  # Decompressed GTF file
    output:
        star="star_output/featurecounts.tsv",  # Read count tsv file (STAR version)
        hisat2="hisat2_output/featurecounts.tsv"  # Read count tsv file (HISAT2 version)
    run:
        count_star=output.star
        count_hisat2=output.hisat2
        bam_star=input.star
        bam_hisat2=input.hisat2
        end=""
        if len(params.ends) == 2:
            end='-p'
        # Additional featureCounts flags
        # -F: format of the annotation file. 'GTF' by default.
        # -g: attribute type. 'gene_id' by default. 
        # -L: set for long-read inputs 
        # -s: strand-specificity. 0 (unstranded & default), 1 (stranded), 2 (reversely stranded)
        shell("set +o pipefail; "
              "featureCounts {end} "    
              "-s0 "
              "-T {THREADS} "
              "-a {params.gtf} "
              "-o {count_star} "
              "{bam_star} &> {count_star}.log && "
              "featureCounts {end} "    
              "-s0 "
              "-T {THREADS} "
              "-a {params.gtf} "
              "-o {count_hisat2} "
              "{bam_hisat2} &> {count_hisat2}.log")


    
    
