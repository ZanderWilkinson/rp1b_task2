# RP1b Task2: Variant Calling Pipeline
This repository contains all the three scripts (task2_part1.py, task2_part2.py, main.py) and files, needed for the RP1b Task 2 assignment (two real Ecoli FASTQ files too large for upload, found in CLIMB mentioned below). Aimed at exploring genome mutation simulation, perfect read simulation and multi-caller pipelines and validation. This README file contains the complete guide to run these scripts, and includes a discussion of results with a worked example.

## Setup 
This section highlights the dependencies, tools and environment required to run the scripts within this project. The workflow of the pipeline is directed by the main Python script, main.py, which chronologically executes functions, and BASH scripts within Python using the subprocess module. 

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
This discussion is structured linearly in the format of main.py, providing examples of important files which are generated as the pipeline runs and addresses insights into the performance / functionality of the pipeline. 

**1. Running Simulation**

As the pipelines runs it generates mutations on both the E.Coli and NC genomes, these mutations (ground truth) are stored in CSV files (Ecoli / NC _Ground_Truth.csv). 100bp reads are then generated on both genomes, which are written to two FASTQ files (Ecoli / NC _mutated_reads.fastq), and these are used in a bcftools variant caller pipelines for single end reads. The precision and recall on both these VCFs are calculated using the equations below. 
   
<p align="center"> 

<img width="204" height="155" alt="image" src="https://github.com/user-attachments/assets/9366ae6e-1776-4d3f-a56c-a70c5ff130eb" /> 
  
<img width="190" height="71" alt="image" src="https://github.com/user-attachments/assets/d62fbd73-5f39-4349-8943-d0d5c3584b01" />

</p>

**2. Precision and Recall**

| Value | Definition | 
| ----------- | ----------- |
| True Positive (TP) | A true mutation that was correctly identified by the caller | 
| False Positive (FP) | A variant called which is not a true mutation | 
| False Negative (FN) | A true mutation that was missed by the caller | 

Precison measures the quality of postive calls =  TP / (TP + FP)

Recall measures the quantity of postive calls = TP / (TP + FN)

From the worked example below we can see that there is a high percentage for both precision and recall on the two genomes. This is due to the use of perfect reads and high coverage, as most true variants are covered. There is some error introduced to the variant callers due to left-alignment, where it reports an indel using an anchor base that is shifted too far to the right. This results in a simple indel being reported as a shifted change (e.g. ATGT -> ATGTGGTG). This causes a discrepancy between the ground truth and VCF, in reality both precision and recall should be 100% if this alignment issue was fixed with hard coding the offsets. However, it is important to truthfully report on issues faced with variant callers. 


<p align="center">  <img width="364" height="69" alt="image" src="https://github.com/user-attachments/assets/dce4379c-255a-489d-89c6-511c0de6a073" /> </p>


**3. Multi-Caller Pipeline**

Firstly, the multi-caller pipeline was ran on our simulated data by creating paired-end reads on the Ecoli genome. The precision and recall for all three VCFs was then scored, as shown below: 

<p align="center"> <img width="346" height="80" alt="image" src="https://github.com/user-attachments/assets/1d722916-f84d-4687-9cc6-70ea2aa0a8d7" />  </p>

The Merged VCF focuses on maximising the number of detected variants (Recall) by combining both VCFs from BCF and Snippy. It recall is consistenly higher than the two individual callers, confirming the merging processes ability to capture all true variants found by either tool. However, during the merge, False Postive calls are inherented from either caller resulting in low precision, in this case the BCF low precision. 

Secondly, the multi-caller pipeline was then run on a real pair of Ecoli reads. Since their is no ground truth, precision and recall cannot be calculated. Therefore, a trust-score system was used to assign confidence ratings and quality scores. An example of this is below, but the full report can be viewed in the CSV file: 

| Trust Score | Meaning | 
| ----------- | ----------- |
| 100 | Consensus between callers, high confidence call | 
| 50 | Either one of the callers, medium confidence call | 

A summary of total variants called is then printed to the terminal, including the number of high/medium scores. From this example we can see the majority of variants are assigned a high confidence score, meaning both callers agreed on the exact variant and position. As previously seen, snippy seems to be much more precise as it uses an optimised calling approach. Due to this, it is likely that the medium confidence calls are false postives inherited from the less precise BCF VCfs. 

<p align="center">  <img width="303" height="68" alt="image" src="https://github.com/user-attachments/assets/25e32f27-459a-41ba-8964-a9c1154e98bb" />  </p>


## Tview Evaluation
To externally validate the trust score, a manual inspection of variants can be done in a genome viewer (tview). Run the following command in the virtual enviroment. To visually confirm a trust score of 100, we should see an agreement of reads showing the variant base, with high quality score and little to no noise. A trust score of 50 may be viewed with a mix of signals in the reads, including low coverage, low quality or disagreeing reads. In the example below we can see that between position 58-59 the alignment switches to a C, and virtually every read, aligned to this position, agree with the call. This visualisation is in agreement with our CSV table which states this variant is a high trust call. Furthermore, in this image around postion 47 their is a missing single base insertion, as well as discrepancies in quality scores, this was indicated in our trust score CSV detecting it as a medium trust variant.

```samtools tview real_bcftools.sorted.bam EcoliK12-MG1655.fasta```

<img width="396" height="1007" alt="image" src="https://github.com/user-attachments/assets/0d6cdf2c-8799-4fb0-b5b1-fc5c9b3a6b6a" />






