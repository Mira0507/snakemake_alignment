
#################################### Defined by users #################################
configfile: "config/config_single.yaml"    # Sets path to the config file
#######################################################################################

shell.prefix('set -euo pipefail; ')
shell.executable('/bin/bash')



rule all: 
    input: 
        expand("fastq/{sample}_{end}.fastq.gz", sample=list(config['SAMPLE'].keys()), end=config['END']), # FASTQ files
        expand("reference/{ref}", ref=config['REFERENCE'][1:]),  # Reference genome and annotation (GTF) files
        expand("hisat2_output/{sample}.bam", sample=list(config['SAMPLE'].keys())),   # HISAT2 output BAM files
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", sample=list(config['SAMPLE'].keys()))   # STAR output BAM files




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
        shell("fastq-dump --split-files {sra} --gzip -X 100000")   # -X is for testing
        for i in range(len(output)):
            i += 1
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
        shell("set +o pipefail; " 
              "wget -c {link}.gz -O {output}.gz && " 
              "gzip -d {output}.gz")


rule index_hisat2:
    """
    This rule constructs HISAT2 index files
    """
    input: 
        expand("reference/{gen}", gen=config['REFERENCE'][1])    # Decompressed reference genome file
    output:
        expand("reference/hisat2_index/hisat2_index.{number}.ht2", number=[x+1 for x in range(8)])   # HISAT2 indexing files
    threads: 8
    shell:
        "set +o pipefail; "
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
        "reference/star_index/Genome",   # STAR indexing files
        "reference/star_index/SA",       # STAR indexing files
        "reference/star_index/SAindex"   # STAR indexing files
    threads: 8
    shell:
        "set +o pipefail; "
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
        index=expand("reference/hisat2_index/hisat2_index.{number}.ht2", number=[x+1 for x in range(8)])  # HISAT2 indexing files
    output:
        temp("hisat2_output/{sample}.sam"),
        "hisat2_output/{sample}.bam"    # Bam files
    params:
        read=config["END"],  # e.g. [1] for single-end, [1, 2] for paired-end
        ext=config['FASTQ_EXT'],       # extension of the FASTQ files (e.g. .fastq.gz)
        indexing=config['INDEX_HISAT'] # HISAT2 indexing location and file name prefix
    threads: 8
    run:
        r1="fastq/" + wildcards.sample + "_1" + params.ext
        r2=""
        read="-U " + r1
        if len(params.read) == 2: 
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
        index1="reference/star_index/Genome",
        index2="reference/star_index/SA",
        index3="reference/star_index/SAindex"
    output:
        "star_output/{sample}Aligned.sortedByCoord.out.bam"     # Bam files
    params:
        read=config["END"],  # e.g. [1] for single-end, [1, 2] for paired-end
        indexing=config["INDEX_STAR"],  # STAR indexing file directory
        ext=config['FASTQ_EXT']         # extension of the FASTQ files (e.g. fastq.gz)
    threads: 8
    run:
        r1= "fastq/" + wildcards.sample + "_1" + params.ext + " " 
        r2=""
        if len(params.read) == 2: 
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
                    "--chimOutType Junctions")     
