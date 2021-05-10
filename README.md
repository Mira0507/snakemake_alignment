### Creating Snakemake Workflow for STAR & HISAT2 Alignments 

#### 1. Conda environment

- [config/conda_env.yml](https://github.com/Mira0507/snakemake_alignment/blob/master/config/conda_env.yml)

#### 2. Snakemake 

- [Snakefile](https://github.com/Mira0507/snakemake_alignment)

- [config/config_single.yaml (single-end testing)](https://github.com/Mira0507/snakemake_alignment/blob/master/config/config_single.yaml)

- [config/config_paired.yaml (paired-end testing)](https://github.com/Mira0507/snakemake_alignment/blob/master/config/config_paired.yaml)

#### 3. Differential expression (DE) analysis

- [DE_analysis/DE.Rmd](): R script

- [config/config_single.R](): R config file (single-end reads)

- [config/config_paired.R](): R config file (paired-end reads)

- [config/sample_single.csv](): sample table csv file (single-end reads)

- [config/sample_paired.csv](): sample table csv file (paired-end reads)

- [config/conda_r.yaml](): conda environment for the rule running R script

#### 4. Error 

- [errors](https://github.com/Mira0507/snakemake_alignment/tree/master/errors): error nots 


#### 5. Running snakemake

- Reference: [Snakemake Command Line Arguments](https://snakemake.readthedocs.io/en/stable/executing/cli.html)

**- Dry run**


```bash

#!/bin/bash

snakemake -n

```


**- DAG visualization**

```bash

#!/bin/bash


# The dot commend requires graphviz (downloadable via conda)
snakemake --dag | dot -Tpdf > dag.pdf

```


**- Run**

```bash
#!/bin/bash

# Either -j or --cores assignes the number of cores
# --use-conda is needed for the rule running DE analysis with isolated conda env
snakemake -j 8 --use-conda

```
