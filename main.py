from task2_part1 import * 
from task2_part2 import *

if __name__ == '__main__':
    
    # PART 1
    #-----------------------------------------------------------------------------------------------------------
    ecoli_genome = get_fasta("EcoliK12-MG1655.fasta")
    nc_genome = get_fasta("NC_037282.1.fasta")
    
    print(f"E. coli genome length: {len(ecoli_genome)} bp")
    print(f"NC genome length: {len(nc_genome)} bp")
    print("---------------------------------------------------------------------------------------------------")
    
    ecoli_indels = indel_index(ecoli_genome)
    nc_indels = indel_index(nc_genome)
    
    print(f"E. coli: {len(ecoli_indels)} mutations")
    print(f"NC: {len(nc_indels)} mutations")
    print("---------------------------------------------------------------------------------------------------")
    
    ecoli_mutated = mutate_indels(ecoli_genome, ecoli_indels)
    nc_mutated = mutate_indels(nc_genome, nc_indels)      
    
    print("Mutations Created")
    print("Writing Mutations To CSV")
    print("---------------------------------------------------------------------------------------------------")   
          
    mutations_csv(ecoli_indels, "EcoliK12-MG1655", "Ecoli_GroundTruth.csv")
    mutations_csv(nc_indels, "NC_037282.1", "NC_GroundTruth.csv")
    
    print("Mutation CSV's Written")
    print("Making Reads")
    print("---------------------------------------------------------------------------------------------------")
          
    ecoli_mutated_reads = simulate_reads(ecoli_mutated)
    nc_mutated_reads = simulate_reads(nc_mutated)
    
    print(f"E. coli: {len(ecoli_mutated_reads)} reads")
    print(f"NC: {len(nc_mutated_reads)} reads")
    print("---------------------------------------------------------------------------------------------------")
    
    print("Writing Reads To FASTQ Files")
    print("---------------------------------------------------------------------------------------------------")
          
    write_reads_to_fastq(ecoli_mutated_reads, "Ecoli_mutated_reads.fastq")
    write_reads_to_fastq(nc_mutated_reads, "Nc_mutated_reads.fastq")
    
    print("Tasks Complete: Generated Mutation CSV's -- Read's Generated -- FASTQ Files Generated ")
    print("---------------------------------------------------------------------------------------------------")
    print("Now Calling Variants on Ecoli and NC")
    print("---------------------------------------------------------------------------------------------------")
    
    ecoli_vcf = call_variants("EcoliK12-MG1655.fasta","Ecoli_mutated_reads.fastq","Ecoli_Variants")    
    nc_vcf = call_variants("NC_037282.1.fasta","Nc_mutated_reads.fastq","Nc_Variants")    
    
    print("VCF Files created:")
    print("Ecoli_Variants.vcf")
    print("Nc_Variants.vcf")    
    print("---------------------------------------------------------------------------------------------------")
          
    print("Scoring Variants and Generating Comparison CSVs")
    print("---------------------------------------------------------------------------------------------------")    
    
    ecoli_metrics = score_variants(ecoli_genome, ecoli_indels, ecoli_vcf,"Ecoli_Variants")
    ecoli_precision = ecoli_metrics['Precision']
    ecoli_recall = ecoli_metrics['Recall']    
   
    nc_metrics = score_variants(nc_genome, nc_indels, nc_vcf, "Nc_Variants")
    nc_precision = nc_metrics['Precision']
    nc_recall = nc_metrics['Recall']
    
    print("Scoring Complete. Comparison CSV files generated.")
    print("---------------------------------------------------------------------------------------------------")
    print("--- Precision and Recall Summary ---")
    print(f"E. coli Precision: {ecoli_precision:.4f} | Recall: {ecoli_recall:.4f}")
    print(f"NC Precision: {nc_precision:.4f} | Recall: {nc_recall:.4f}")
    
    print("---------------------------------------------------------------------------------------------------")
    print("Part 1 Complete") 
    print("---------------------------------------------------------------------------------------------------")
    
    # PART 2
    #--------------------------------------------------------------------------------------------------------------
    
    print("Starting Part 2")
    print("---------------------------------------------------------------------------------------------------")
    print("Making Paired End Reads On Mutated Ecoli Genome")
    print("---------------------------------------------------------------------------------------------------")  
          
    ecoli_reverse_mutated_reads = reverse_complement_reads(ecoli_mutated_reads)
          
    print("Reads Made")
    print("---------------------------------------------------------------------------------------------------")
    print("Writing to FASTQ File")
    print("---------------------------------------------------------------------------------------------------")
          
    write_reads_to_fastq(ecoli_reverse_mutated_reads, "Ecoli_reverse_mutated_reads.fastq")
          
    print("FASTQ Generated")
    print("---------------------------------------------------------------------------------------------------")
          
    ref = "EcoliK12-MG1655.fasta"
    
    read1_sim = "Ecoli_mutated_reads.fastq"
    read2_sim = "Ecoli_reverse_mutated_reads.fastq"  
    
    read1_real = "SRR25083113_1.fastq.gz"  
    read2_real = "SRR25083113_2.fastq.gz" 
    
    print("Running BCF and Snippy on Mutated Ecoli Paired End Reads")
    print("Merging The Two VCFs")
    print("---------------------------------------------------------------------------------------------------")
          
    results_sim = run_both_callers(ref, read1_sim, read2_sim, "simulation")   
    
    print("VCFs Generated")
    print("---------------------------------------------------------------------------------------------------")
    print("Scoring Precion and Recall for BCF, Snippy and Merged VCFs")
    print("---------------------------------------------------------------------------------------------------")
    
    ecoli_vcf_bcf = results_sim['bcftools_vcf']
    ecoli_vcf_snippy = results_sim['snippy_vcf']
    ecoli_vcf_merged = results_sim['merged_vcf'] 

    ecoli_metrics_bcf = score_variants(ecoli_genome, ecoli_indels, ecoli_vcf_bcf,"Ecoli_Variants_BCF")
    ecoli_metrics_snippy = score_variants(ecoli_genome, ecoli_indels, ecoli_vcf_snippy,"Ecoli_Variants_Snippy")
    ecoli_metrics_merged = score_variants(ecoli_genome, ecoli_indels, ecoli_vcf_merged,"Ecoli_Variants_Merged")
    
    ecoli_precision_bcf = ecoli_metrics_bcf['Precision']
    ecoli_recall_bcf = ecoli_metrics_bcf['Recall']   
    
    ecoli_precision_snippy = ecoli_metrics_snippy['Precision']
    ecoli_recall_snippy = ecoli_metrics_snippy['Recall']   
    
    ecoli_precision_merged = ecoli_metrics_merged['Precision']
    ecoli_recall_merged = ecoli_metrics_merged['Recall']   
  
    print("--- Precision and Recall Summary ---")
    print(f"BCF Precision: {ecoli_precision_bcf:.4f} | Recall: {ecoli_recall_bcf:.4f}")
    print(f"Snippy Precision: {ecoli_precision_snippy:.4f} | Recall: {ecoli_recall_snippy:.4f}")
    print(f"Merged Precision: {ecoli_precision_merged:.4f} | Recall: {ecoli_recall_merged:.4f}")
    print("---------------------------------------------------------------------------------------------------")
    
    print("Running BCF, Snippy On Real Ecoli Paired End Reads")
    print("Merging The Two VCFs")
    print("---------------------------------------------------------------------------------------------------")
    
    results_real = run_both_callers(ref, read1_real, read2_real, "real")
    
    print("VCFS Made")
    print("Generate Scoring For The VCFS")
    print("---------------------------------------------------------------------------------------------------")
    
    trust_metrics_real = annotate_and_report_trust(results_real, "Ecoli_Real")
    
    print("Scoring Complete")
    print(f"Report written to Ecoli_Real_Trust_Report.csv")
    print("---------------------------------------------------------------------------------------------------")
    print(f"Total Variants Called: {trust_metrics_real['Total']}")
    print(f"High Confidence (TRUST=100): {trust_metrics_real.get('100', 0)}")
    print(f"Medium Confidence (TRUST=50): {trust_metrics_real.get('50', 0)}")
    print("---------------------------------------------------------------------------------------------------")
    print("Part 2 Complete")
   


