### Creating Snakemake Workflow for STAR & HISAT2 Alignments 

#### 1. Conda environment

- [config/conda_env.yml](https://github.com/Mira0507/snakemake_alignment/blob/master/config/conda_env.yml)

#### 2. Snakemake 

- [Snakefile](https://github.com/Mira0507/snakemake_alignment): snakefile 

- [config/config_single.yaml](https://github.com/Mira0507/snakemake_alignment/blob/master/config/config_single.yaml): snakemake config file (single-end testing)

- [config/config_paired.yaml](https://github.com/Mira0507/snakemake_alignment/blob/master/config/config_paired.yaml): snakemake config file (paired-end testing)

#### 3. Differential expression (DE) analysis

- [DE_analysis/DE.Rmd](https://github.com/Mira0507/snakemake_alignment/blob/master/DE_analysis/DE.Rmd): R script

- [config/config_single.R](https://github.com/Mira0507/snakemake_alignment/blob/master/config/config_single.R): R config file (single-end testing)

- [config/config_paired.R](https://github.com/Mira0507/snakemake_alignment/blob/master/config/config_paired.R): R config file (paired-end testing)

- [config/sample_single.csv](https://github.com/Mira0507/snakemake_alignment/blob/master/config/sample_single.csv): sample table csv file (single-end testing)

- [config/sample_paired.csv](https://github.com/Mira0507/snakemake_alignment/blob/master/config/sample_paired.csv): sample table csv file (paired-end testing)

- [config/conda_r.yaml](https://github.com/Mira0507/snakemake_alignment/blob/master/config/conda_r.yaml): conda environment for the rule running R script

#### 4. Error 

- [errors](https://github.com/Mira0507/snakemake_alignment/tree/master/errors): error nots 


#### 5. Running snakemake

- Reference: [Snakemake Command Line Arguments](https://snakemake.readthedocs.io/en/stable/executing/cli.html) (paired-end testing)

- Modifying files:
    - **Snakefile**: configfile: "path/to/config.yaml" (path to the your snakemake config file)
    - **config/config.yaml**: snakemake config file
    - **DE_analysis/DE.Rmd**: source("path/to/config.R") 
    - **config/config.R**: species, sample.csv (file path), alpha (max FDR), mLog (minimum log2FoldChange of interest)
    - **config/sample.csv**: the sample column has to correspond to the sample info in your snakemake config file

- **Dry run**


```bash

#!/bin/bash

snakemake -n

```


- **DAG visualization**

```bash

#!/bin/bash


# The dot commend requires graphviz (downloadable via conda)
snakemake --dag | dot -Tpdf > dag.pdf

```


- **Run**

```bash
#!/bin/bash

# Either -j or --cores assignes the number of cores
# --use-conda is needed for the rule running DE analysis with isolated conda env
snakemake -j 10 --use-conda

```
