
#################################### Defined by users #################################
configfile: "config/config_paired1.yaml"    # Sets path to the config file
#######################################################################################



rule all: 
    input: 
        expand("reference/{ref}", ref=config['REFERENCE'][1:]),  # Reference genome and annotation (GTF) files
        expand("hisat2_output/{sample}.bam", sample=list(config['SAMPLE'].keys())),   # Output BAM files (HISAT2)
        expand("star_output/{sample}.bam", sample=list(config['SAMPLE'].keys())),     # Output BAM files (STAR)
        expand("{aligner}_output/featurecounts.tsv", aligner=config['ALIGNER']), # featureCount output (read count matrix) 
        expand("DE_analysis/{outfile}", outfile=config['DE_OUTPUT'])             # DE analysis result



rule get_fastq:   
    """
    This rule downloads SRA and converts to FASTQ files
    """
    output:
        expand("fastq/{{sample}}_{end}.fastq.gz", end=config['END'])
    params:
        dic=config['SAMPLE']
    run:
        sra=params.dic[wildcards.sample]
        shell("fastq-dump --split-files {sra} --gzip")   # -X can be added for testing
        for i in range(1, len(output)+1):
            shell("mv {sra}_{i}.fastq.gz fastq/{wildcards.sample}_{i}.fastq.gz")


rule get_reference:    
    """
    This rule downloads and decompresses reference files
    """
    params:
        reflink=config['REFERENCE'][0]
    output:
        "reference/{ref}"  # Decompressed reference files
    run:
        link=params.reflink + wildcards.ref
        shell("wget -c {link}.gz -O {output}.gz && " 
              "gzip -d {output}.gz")


rule index_hisat2:
    """
    This rule constructs HISAT2 index files
    """
    input: 
        expand("reference/{gen}", gen=config['REFERENCE'][1])    # Decompressed reference genome file
    output:
        expand("{path}.{number}.ht2", path=config['INDEX_HISAT'], number=[x for x in range(1, 9)])   # HISAT2 indexing files
    threads: 16
    shell:
        "hisat2-build -f -o 4 "
        "-p {threads} "
        "--seed 67 "
        "{input} "
        "hisat2_index && "
        "mv *.ht2 reference/hisat2_index"


rule index_star:
    """
    This rule constructs STAR index files
    """
    input:
        fa=expand("reference/{gen}", gen=config['REFERENCE'][1]),  # Decompressed reference genome file
        gtf=expand("reference/{anno}", anno=config['REFERENCE'][2])  # Decompressed GTF file
    output:
        expand("reference/star_index/{index}", index=config['INDEX_STAR'][1:])
    threads: 16
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir reference/star_index "
        "--genomeFastaFiles {input.fa} "
        "--sjdbGTFfile {input.gtf}"


rule align_hisat2:    # Creates bam files in hisat2_output directory"
    """
    This rule aligns the reads using HISAT2    
    """
    input:
        fastq=expand("fastq/{{sample}}_{end}.fastq.gz", end=config["END"]),   # Gzipped FASTQ files
        index=expand("{path}.{number}.ht2", path=config['INDEX_HISAT'], number=[x for x in range(1, 9)])  # HISAT2 indexing files
    output:
        temp("hisat2_output/{sample}.sam"),
        "hisat2_output/{sample}.bam"    # Bam files
    params:
        ext=config['FASTQ_EXT'],       # extension of the FASTQ files (e.g. .fastq.gz)
        indexing=config['INDEX_HISAT'] # HISAT2 indexing location and file name prefix
    threads: 16
    run:
        r1="fastq/" + wildcards.sample + "_1" + params.ext
        r2=""
        read="-U " + r1
        if len(input.fastq) == 2:    # if paired-end
            r2= "fastq/" + wildcards.sample + "_2" + params.ext  
            read="-1 " + r1 + " -2 " + r2
        shell("hisat2 -q -p {threads} "
              "--seed 23 "
              "--summary-file hisat2_output/summary_{wildcards.sample}.txt "
              "-x {params.indexing} "
              "{read} "
              "-S {output[0]} && "
              "samtools view -bS "
              "-@ {threads} "
              "{output[0]} > {output[1]}")



rule align_star:   # Creates bam files in star_output directory"
    """
    This rule aligns the reads using STAR two-pass mode
    """
    input:
        gtf=expand("reference/{anno}", anno=config['REFERENCE'][2]),   # Decompressed GTF file
        fastq=expand("fastq/{{sample}}_{end}.fastq.gz", end=config['END']),    # Gzipped FASTQ files
        index=expand("reference/star_index/{index}", index=config['INDEX_STAR'][1:]) # STAR indexing files
    output:
        "star_output/{sample}.bam"     # Bam files
    params:
        indexing=config["INDEX_STAR"][0],  # STAR indexing file directory
        ext=config['FASTQ_EXT']         # extension of the FASTQ files (e.g. fastq.gz)
    threads: 16
    run:
        r1= "fastq/" + wildcards.sample + "_1" + params.ext + " " 
        r2=""
        if len(input.fastq) == 2:   # if paired-end
            r2= "fastq/" + wildcards.sample + "_2" + params.ext + " " 
        shell("STAR --runThreadN {threads} "  
                "--runMode alignReads "  
                "--readFilesCommand zcat "
                "--genomeDir {params.indexing} " 
                "--sjdbGTFfile {input.gtf} "  
                "--sjdbOverhang 100 "  
                "--readFilesIn {r1}{r2}"  
                "--outFileNamePrefix star_output/{wildcards.sample} "
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
                "--chimOutType Junctions && "
              "mv star_output/{wildcards.sample}Aligned.sortedByCoord.out.bam {output}")   


rule featurecounts:
    """
    This rule assesses read counts using featureCounts
    """
    input:
        bam=expand("{{aligner}}_output/{sample}.bam", sample=list(config['SAMPLE'].keys())),  # BAM files 
        gtf=expand("reference/{anno}", anno=config['REFERENCE'][2])  # Decompressed GTF file
    params:
        ends=config['END']   # Read ends (e.g. [1] or [1, 2])
    output:
        "{aligner}_output/featurecounts.tsv"  # Read count tsv file 
    threads: 8
    run:
        r=""
        if len(params.ends) == 2:   # if paired-end
            r='-p'
        # Additional featureCounts flags
        # -F: format of the annotation file. 'GTF' by default.
        # -g: attribute type. 'gene_id' by default. 
        # -L: set for long-read inputs 
        # -s: strand-specificity. 0 (unstranded & default), 1 (stranded), 2 (reversely stranded)
        shell("featureCounts {r} "    
              "-s0 "
              "-T {threads} "
              "-a {input.gtf} "
              "-o {output} "
              "{input.bam} &> {output}.log")


rule deseq2:
    """
    This rule performs differential expression (DE) analysis using DESeq2 in R
    """
    input:
        reads=expand("{aligner}_output/featurecounts.tsv", aligner=config["ALIGNER"]),  # Assessed read counts 
        rmd="DE_analysis/DE.Rmd",                    # R scripts for downstream DE analysis
        rconfig=config["RCONFIG"]                    # config file for Rmd 
    output:
        expand("DE_analysis/{outfile}", outfile=config['DE_OUTPUT']),     # Output files from Rmd 
        temp("DE_analysis/trf_input.fa"),
        temp("DE_analysis/DE.md")
    conda:
        config['CONDA']     # conda env config file
    shell:
        "Rscript -e \"rmarkdown::render('{input.rmd}')\""   # Requires "snakemake -j 10 --use-conda" command when running snakemake
