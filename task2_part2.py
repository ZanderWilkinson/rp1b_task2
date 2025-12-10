# Import modules
import shutil 
import csv
import subprocess
import os
from collections import defaultdict


# Call variants using minimap2 and bcftools with paired-end reads
def call_variants_bcftools(ref, read1, read2, out_prefix):
    
    # Define output file paths
    sorted_bam = f"{out_prefix}.sorted.bam"
    vcf = f"{out_prefix}.vcf" 
    
    # Align paired-end reads to reference using minimap2
    minimap_proc = subprocess.Popen(["minimap2", "-ax", "sr", ref, read1, read2], stdout=subprocess.PIPE)    
    # Convert SAM to BAM format
    view_proc = subprocess.Popen(["samtools", "view", "-b"], stdin=minimap_proc.stdout, stdout=subprocess.PIPE)
    
    # Sort alignments and write to BAM file
    with open(sorted_bam, "wb") as bam_out:
        subprocess.run(["samtools", "sort", "-o", sorted_bam], stdin=view_proc.stdout, check=True)
    
    # Index the sorted BAM file
    subprocess.run(["samtools", "index", sorted_bam], check=True)
    
    # Generate pileup and call variants with bcftools
    with open(vcf, "w") as vcf_out:
        mpileup_proc = subprocess.Popen(["bcftools", "mpileup", "-f", ref, sorted_bam], stdout=subprocess.PIPE)
        subprocess.run(["bcftools", "call", "-mv", "-Ov"], stdin=mpileup_proc.stdout, stdout=vcf_out, check=True)
    
    return vcf, sorted_bam


# Call variants using snippy with paired-end reads
def call_variants_snippy(ref, read1, read2, out_prefix):
    
    # Define snippy output directory
    output_dir = f"{out_prefix}_snippy_dir"
    # Build snippy command with all parameters
    cmd = ["snippy","--cpus", "2","--outdir", output_dir,"--ref", ref,"--R1", read1,"--R2", read2,"--force" ]
    
    # Execute snippy and handle errors
    try:
        result = subprocess.run(cmd, shell=False, check=True, capture_output=True, text=True)
        # Snippy outputs VCF as snps.vcf in the output directory
        vcf = os.path.join(output_dir, "snps.vcf")
        
        # Check if VCF was successfully created
        if os.path.exists(vcf):            
            return vcf, output_dir
        else:            
            return None, output_dir
    
    # Handle subprocess errors        
    except subprocess.CalledProcessError as e:
        print("Error: Snippy failed to run.")
        raise

# Merge two VCF files using bcftools merge
def merge_vcfs(vcf1, vcf2, output_vcf):     
    
    # Define compressed VCF file paths
    vcf1_gz = f"{vcf1}.gz"
    vcf2_gz = f"{vcf2}.gz"    
    
    # Compress first VCF with bgzip
    subprocess.run(["bgzip", "-f", "-c", vcf1], stdout=open(vcf1_gz, "wb"), check=True)     
    # Index first compressed VCF
    subprocess.run(["tabix", "-f", "-p", "vcf", vcf1_gz], check=True)     
    # Compress second VCF with bgzip
    subprocess.run(["bgzip", "-f", "-c", vcf2], stdout=open(vcf2_gz, "wb"), check=True)  
    # Index second compressed VCF
    subprocess.run(["tabix", "-f", "-p", "vcf", vcf2_gz], check=True)    
   
    # Merge both VCFs into single output file
    with open(output_vcf, "w") as out:
        subprocess.run(["bcftools", "merge", vcf1_gz, vcf2_gz, "-o", output_vcf, "-O", "v"], check=True)    
    
    return output_vcf        
        
        
# Run both variant callers and merge results
def run_both_callers(ref, read1, read2, out_prefix):
    
    # Run bcftools variant calling pipeline
    bcftools_vcf, bcftools_bam = call_variants_bcftools(ref, read1, read2, f"{out_prefix}_bcftools")     
    # Run snippy variant calling pipeline
    snippy_vcf, snippy_dir = call_variants_snippy(ref, read1, read2, f"{out_prefix}_snippy")    
   
    # Merge VCF files from both callers
    merged_vcf = f"{out_prefix}_merged.vcf"
    merge_vcfs(bcftools_vcf, snippy_vcf, merged_vcf)       
    
    # Return dictionary with all output file paths
    return {
        'bcftools_vcf': bcftools_vcf,
        'bcftools_bam': bcftools_bam,
        'snippy_vcf': snippy_vcf,
        'snippy_dir': snippy_dir,
        'merged_vcf': merged_vcf}


# Generate reverse complement of DNA sequence
def make_reverse_complement(sequence):   
    
    # Define complementary base pairs
    COMPLEMENT_MAP = str.maketrans("ACGT", "TGCA")
    # Get complement of sequence
    complement = sequence.translate(COMPLEMENT_MAP)
    
    # Reverse the complement to get reverse complement
    return complement[::-1]

# Create reverse complement reads for paired-end simulation
def reverse_complement_reads(original_reads): 
    
    # List to store reverse complement reads
    rc_reads = []    
    
    # Process each read
    for read_id, sequence in original_reads:
        
        # Generate reverse complement of sequence
        rc_sequence = make_reverse_complement(sequence)        
        
        # Keep original read ID
        read_id = f"{read_id}"      
       
        # Add reverse complement read to list
        rc_reads.append((read_id, rc_sequence))
        
    return rc_reads

# Parse VCF file and extract variants as a set for comparison
def parse_vcf_to_set(vcf_path):
    
    # Set to store unique variants as (position, ref, alt) tuples
    variant_set = set()
    
    # Open and parse VCF file
    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Split line into columns
                cols = line.strip().split('\t')
                
                # Extract variant information
                try:
                    pos = int(cols[1])
                    ref = cols[3]
                    alts = cols[4].split(',')
                    
                    # Add each alternate allele as separate variant
                    for alt in alts:
                        variant_set.add((pos, ref, alt))
                
                # Skip malformed lines        
                except (IndexError, ValueError):
                    continue
    
    # Return empty set if file not found            
    except FileNotFoundError:
        return set()

    return variant_set

# Annotate merged VCF with trust scores based on caller agreement
def annotate_trust_score(bcftools_vcf, snippy_vcf, merged_vcf):
    
    # Parse variants from each caller into sets
    bcftools_set = parse_vcf_to_set(bcftools_vcf)
    snippy_set = parse_vcf_to_set(snippy_vcf)
    
    # Create temporary file for annotated VCF
    temp_vcf = merged_vcf + ".temp_annotated"
    
    # Define new INFO field headers
    trust_info = '##INFO=<ID=TRUST,Number=1,Type=Integer,Description="Consensus Score (100=both, 50=one caller)">\n'
    callers_info = '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Callers reporting this variant">\n'

    # Read input VCF and write annotated output
    with open(merged_vcf, 'r') as vcf_in, open(temp_vcf, 'w') as vcf_out:
        
        # Track if new INFO headers have been written
        header_written = False
        
        # Process each line of VCF
        for line in vcf_in:
            # Write existing header lines
            if line.startswith('##'):
                vcf_out.write(line)
                
                # Insert new INFO headers after bcftools_merge line
                if line.startswith('##bcftools_merge') and not header_written:
                    vcf_out.write(trust_info)
                    vcf_out.write(callers_info)
                    header_written = True
                continue

            # Write column header line and ensure INFO headers are added
            if line.startswith('#CHROM'):
                if not header_written: 
                    vcf_out.write(trust_info)
                    vcf_out.write(callers_info)
                vcf_out.write(line)
                continue
            
            
            # Parse variant record
            cols = line.strip().split('\t')
            
            # Extract position, ref, alt, and INFO fields
            try:
                pos, ref, alts_str, info = int(cols[1]), cols[3], cols[4], cols[7]
            except (IndexError, ValueError):
                vcf_out.write(line)
                continue
            
            # Process each alternate allele
            alts = alts_str.split(',')
            trust_score = 0
            callers = []
            
            # Check which callers reported each variant
            for alt in alts:
                
                # Create variant tuple for lookup
                variant_tuple = (pos, ref, alt)
                is_bcftools = variant_tuple in bcftools_set
                is_snippy = variant_tuple in snippy_set
                
                # Track which callers found this variant
                if is_bcftools: callers.append('BCF')
                if is_snippy: callers.append('SNIPPY')
                
                # Assign trust score based on caller agreement
                if is_bcftools and is_snippy:
                    trust_score = 100
                    break 
                elif trust_score < 50 and (is_bcftools or is_snippy):
                    trust_score = 50

            # Add TRUST and CALLERS annotations to INFO field
            if trust_score > 0:
                unique_callers = ','.join(sorted(list(set(callers))))
                info += f";TRUST={trust_score};CALLERS={unique_callers}"
            
            # Write annotated record
            cols[7] = info
            vcf_out.write('\t'.join(cols) + '\n')
    
    # Replace original file with annotated version        
    shutil.move(temp_vcf, merged_vcf)

    return merged_vcf

# Extract trust scores from annotated VCF and write to CSV report
def write_trust_scores_to_csv(vcf_path, output_filename):
    
    # Dictionary to count variants by trust score
    trust_scores = {'100': 0, '50': 0, 'Total': 0}
    # Define CSV columns
    fieldnames = ["POS", "REF", "ALT", "QUAL", "FILTER", "TRUST_SCORE", "CALLERS"]

    # Open VCF and CSV files
    try:
        with open(vcf_path, 'r') as vcf_in, \
             open(output_filename, 'w', newline='') as csv_out:
            
            # Initialize CSV writer
            writer = csv.DictWriter(csv_out, fieldnames=fieldnames)
            writer.writeheader()
            
            # Process each variant in VCF
            for line in vcf_in:
                # Skip header lines
                if line.startswith("#"):
                    continue

                # Count total variants
                trust_scores['Total'] += 1
                
                # Parse variant fields
                cols = line.strip().split("\t")
                pos, ref, alt, qual, filter_status, info = cols[1], cols[3], cols[4], cols[5], cols[6], cols[7]
                
                # Initialize trust score and callers as not available
                trust_score = 'N/A'
                callers_list = 'N/A'
                
                # Parse INFO field into dictionary
                info_fields = dict(item.split("=") for item in info.split(";") if "=" in item)
                
                # Extract TRUST score if present
                if 'TRUST' in info_fields:
                    trust_score = info_fields['TRUST']
                    # Count variants by trust score
                    if trust_score in trust_scores:
                        trust_scores[trust_score] += 1
                
                # Extract CALLERS list if present        
                if 'CALLERS' in info_fields:
                    callers_list = info_fields['CALLERS'].replace(',', ';')
                
                # Write variant information to CSV
                writer.writerow({
                    "POS": pos, "REF": ref, "ALT": alt, "QUAL": qual, 
                    "FILTER": filter_status, "TRUST_SCORE": trust_score, 
                    "CALLERS": callers_list
                })

    # Return empty counts if file not found
    except FileNotFoundError:
        return trust_scores
        
    return trust_scores

# Annotate merged VCF with trust scores and generate CSV report
def annotate_and_report_trust(run_results, out_prefix):
   
    # Extract VCF paths from run results
    bcftools_vcf = run_results['bcftools_vcf']
    snippy_vcf = run_results['snippy_vcf']
    merged_vcf = run_results['merged_vcf'] 
    # Annotate merged VCF with trust scores
    annotated_vcf_path = annotate_trust_score(bcftools_vcf, snippy_vcf, merged_vcf)
    # Generate CSV report with trust scores
    output_filename = f"{out_prefix}_Trust_Report.csv"
    trust_metrics = write_trust_scores_to_csv(annotated_vcf_path, output_filename) 
    
    return trust_metrics