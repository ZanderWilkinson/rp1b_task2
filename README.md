# RP1b Task2: Variant Calling Pipeline
This repository contains all the three scripts (task2_part1.py, task2_part2.py, main.py) and files, needed for the RP1b Task 2 assignment. Aimed at exploring genome mutation simulation, perfect read simulation and multi-caller pipelines and validation. This README file contains the complete guide to run these scripts, and includes a discussion of results with a worked example.

## Setup 
This section highlights the dependencies, tools and enviroment required to run the scripts within this project. The workflow of the pipeline is directed by the main Python script, main.py, which chronologically executes functions, and BASH scripts within Python using the subprocess module. 

### Conda Environment 
This codebase was run in a Conda virutal enviroment within CLIMB computing (path: /shared-team/people/zander/Task2/), due to issues with compatability between bioinformatics tools, specific versions are needed to run these scripts. In the BASH terminal in CLIMB, use the following code to create and activate the necessary environment:

```conda create -n task_2 snippy samtools=1.9```

```conda activate task_2```


### Requirements 
The following modules should be installed within your task_2 environment: 

- Python 3
- random, subprocess, csv, os , shutil, collections (defaultdict)

The following input files should be present in you working directory: 


| Filename | Description | 
| ----------- | ----------- |
| EcoliK12-MG1655.fasta| Ecoli Reference Genome for Mutations and Alignment | 
| NC_037282.1.fasta | Plasmodium Reference Genome | 
| SRR25083113_1.fastq.gz | Real Illumina Paired-end Read 1 | 
| SRR25083113_2.fastq.gz | Real Illumina Paired-end Read 1  | 


### Tools 
This project relies on a number of command-line bioinformatics tools, highlighted in the following table. All tools must be installed within the Conda environment and accesible in your systems working path (they should be automatically installed when creating the environment, but in the case of an error, installs should be checked).


| Tool Name | Description | Purpose |
| ----------- | ----------- |-----------|
| minimap2 | Read alignment | Maps reads to reference genome |
| samtools | Bam Manipulation | Converts SAM to BAM, sorts and indexes alignments |
| bcftools | Variant Calling + Merging | Primary variant caller |
| snippy | Variant Calling | Altnerative variant caller for analysis |
| bgzip | VCF Compression | Required for bcftools merge |
| tabix | VCF Indexing | Required for bcftools merge |

### Running the Script 
The entire analysis can be run by executing the main file once, ensure all required files and scripts are in the same directory and the script is executed in the conda environment:

``` python main.py```


## Part 1: Simulation and Validation 
The functions for part 1 can be viewed in task2_part1.py. They perform mutations, read simulation, variant calling and validation against the ground truth, in the form of precision and recall. The main functions are summarise below. 


| Key Functions | Description | 
| ----------- | ----------- |
| indel_index | Generates a ground truth set of 300 random SNPS and 20 small (1-10bp) indes | 
| mutate_indels | Applies the mutations to the reference genome string | 
| simulate_reads | Generates ~ 30x depth reads of perfect 100bp single end reads from the mutated genome | 
| call_variants | Runs the full minimap2, samtools sort, bcftools mpileup and bcftools call on single reads |
| compare_variants | Compares call variants from VCF against the true mutations | 
| calculate_metrics | Calculates precision and recall from total variants, false postives and true postives |



## Part 2: Multi-Caller Pipeline and Trust Scoring 
The functions for part 2 can be viewed in task2_part2.py. They implement a multi-caller pipeline on both simulated paired end reads and real data, providing a scoring system for the trust of each vcf.

| Key Functions | Description | 
| ----------- | ----------- |
| call_variants_bcftools | Runs the bcftools pipeline for paired-end reads | 
| call_variants_snippy | Runs snippy for variant calling for paired-end reads | 
| merge_vcfs | Uses bgzip, tabix and bcftools merge to combine VCFs from both callers | 
| reverse_complement_reads | Generates a reverse strand reads from simulated mutations, for calling  | 
| annotate_trust_score | Compares variants in merged vcf against bcf and snippy VCFs, giving trust scores | 
| write_trust_scores_to_csv | Generates a final CSV report of trust scores | 



## Discussion + Worked Example 
This discussion is structured linearly in the format of main.py, providing examples of important files which are generated as the pipeline runs and addresses insights into the performance / fuctionality of the pipeline. 

1. Precision and Recall of Simulated Reads



