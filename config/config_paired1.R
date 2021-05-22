my.anno.db <- "EnsDb"  # DB for AnnotationHub 
my.species <- "Homo sapiens"   # Species of the input dataset


featurecounts.path <- list(HISAT2="../hisat2_output/featurecounts.tsv",
                           STAR="../star_output/featurecounts.tsv")  # Paths to featurecounts count matrices 



sample.csv <- "../config/sample_paired1.csv"   # Path to the sample table csv file

alpha <- 0.1   # Set threshold FDR



mLog <- c(-1, 1)  # Set minimum log2 fold change of interest 

basemean.comparison.ylog <- T
lfc.comparison.ylog <- F
padj.comparison.ylog <- F



geneid.species <- "ENSG"  # Ensembl code for gene id given by species (e.g. ENSG (human), ENSMUSG (mouse))


num.dis.genes <- 100      # Number of discordant genes (in percent difference). If num.dis.genes = 100, 100 top- and 100-bottom genes are considered at discordance level between HISAT2 and STAR
