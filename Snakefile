



# This workflow is designed to download fastq files from SRA database. 
# It's possible to perform manually as well (see https://github.com/Mira0507/using_SRA)

#################################### Defined by users #################################
configfile:"config/config.yaml"    # Sets path to the config file

#######################################################################################

THREADS=config["THREADS"]

shell.prefix('set -euo pipefail')
shell.executable('/bin/bash')


rule all:
    input:
        expand("star_output/{sample}Aligned.sortedByCord.out.bam", 
               sample=config["FILE_NAME"])



rule create_dir: 
    """
    This rule creates output directories
    """
    output:
        temp("temp.txt")
    shell:
        "set +o pipefail; "
        "mkdir star_output hisat2_output && touch temp.txt"

rule align_star: 
    """
    This rule aligns the reads using STAR two-pass mode
    """
    input:
        gtf=config["GTF"],
        reads=config["INPUT_PATH"],
        tempfile="temp.txt"
    output:
        "star_output/{sample}Aligned.sortedByCord.out.bam"   
    params:
        indexing=config["INDEX_STAR"],
        files=config["FILE_NAME"]
    run:
        for i in range(len(input.reads)):
            r=input.reads[i]
            p=params.files[i]
            shell("set +o pipefail; "
                  "cd star_output && " 
                  "STAR --runThreadN {THREADS} "  
                    "--runMode alignReads "  
                    "--readFilesCommand zcat "
                    "--genomeDir {params.indexing} " 
                    "--sjdbGTFfile ../{input.gtf} "  
                    "--sjdbOverhang 100 "  
                    "--readFilesIn {r} "  
                    "--outFileNamePrefix {p} "
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
                    "--chimOutType Junctions >> star.log && cd ..")


