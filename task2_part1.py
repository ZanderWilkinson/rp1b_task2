# Import modules
import random
import subprocess
import csv
from collections import defaultdict

# Read FASTA file and return genome sequence as a string
def get_fasta(fasta):
    
    # Initialize list to store genome lines
    genome_list = []
    
    # Open and read the FASTA file
    with open(fasta) as genome_file:
        genome_lines = genome_file.readlines()
        # Skip header lines (starting with '>') and collect sequence lines
        for line in genome_lines:
            if not line.startswith(">"):
                genome_list.append(line)
    # Join all sequence lines and remove newline characters
    genome_string = "".join(genome_list).replace("\n", "")
    
    return genome_string


# Generate random mutations (SNPs and small indels) across the genome
def indel_index(genome_string):
    
    # Define the four DNA bases
    bases = ["A", "C", "G", "T"]
    indels = []
    # Randomly decide how many insertions and deletions
    num_insertions = random.choice(range(1, 20))
    num_deletions = 20 - num_insertions
    # Randomly select positions for insertions
    insertion_indexs = sorted(random.sample(range(len(genome_string)), num_insertions))
    deletion_indexs = []
    snp_indexs = []
    
    # Generate deletion positions that don't overlap with insertions
    while len(deletion_indexs) < num_deletions:
        index = random.randrange(len(genome_string) - 10)
        if index not in insertion_indexs:
            deletion_indexs.append(index)
    
    # Generate 300 SNP positions that don't overlap with indels        
    while len(snp_indexs) < 300:
        index = random.randrange(len(genome_string))
        if index not in deletion_indexs and index not in insertion_indexs and index not in snp_indexs:
            snp_indexs.append(index)
    
    # Create insertion mutations with random lengths (1-9 bp) and sequences        
    for index in insertion_indexs:
        bp_length = random.choice(range(1,10))
        insertion = random.choices(bases, k = bp_length)
        indels.append(("INS", index, bp_length, "".join(insertion)))
    
    # Create deletion mutations with random lengths (1-9 bp)    
    for index in deletion_indexs:
        bp_length = random.choice(range(1,10))
        indels.append(("DEL", index, bp_length, genome_string[index:index+bp_length]))
    
    # Create SNP mutations (single base substitutions)    
    for index in snp_indexs:
        snp = random.choice([base for base in bases if base != genome_string[index]])
        indels.append(("SNP", index, genome_string[index], snp))
    
    # Sort mutations in reverse order by position (highest to lowest)    
    indels.sort(key = lambda x: x[1], reverse=True)    
    
    return indels


# Apply mutations to the genome sequence
def mutate_indels(genome_string,indels):
    
    # Convert genome string to list for easier manipulation
    genome_list = list(genome_string)
    
    # Apply each mutation (processing from end to start preserves earlier positions)
    for type, pos, identifier, change in indels:
        # Insert new bases at the position
        if type == "INS":
            genome_list[pos:pos] = list(change)
        # Replace the base at the position
        elif type == "SNP":
            genome_list[pos] = change
        # Delete bases starting at the position
        elif type == "DEL":
            del genome_list[pos: pos + identifier]
    
    # Convert list back to string        
    return "".join(genome_list)


# Write true mutations to CSV file for later comparison
def mutations_csv(indels, genome_name, filename="Truth_mutations.csv"):
    
    # Sort mutations by position (lowest to highest) for easier reading
    sorted_indels = sorted(indels, key=lambda x: x[1])     
    header = ["NAME", "POS", "REF", "ALT"]

    # Write mutations to CSV file
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        
        # Format each mutation for CSV output
        for name, pos, identifier, change in sorted_indels:
            # Convert to 1-based indexing (VCF standard)
            pos_1 = pos + 1         
            if name == "SNP":
                ref = identifier  # The single REF base
                alt = change      # The single ALT base
            elif name == "DEL":
                ref = f"Length = {identifier}" # The length of the deletion
                alt = f"{change}"   # The sequence that was deleted
            elif name == "INS":
                ref = f"Length = {identifier}" # The length of the insertion
                alt = f"{change}"   # The sequence that was inserted
            else:
                ref = "---"
                alt = "---"

            writer.writerow([name, pos_1, ref, alt])
             
                
# Simulate 30x coverage of 100bp reads from the genome
def simulate_reads(genome_string):
    
    # Calculate number of reads needed for 30x coverage
    number_of_reads = int((30*len(genome_string))/100)
    reads = []
    
    # Generate reads at random positions
    for i in range(number_of_reads):
        start = random.randrange(0, len(genome_string) - 100)
        read = genome_string[start:start+100]
        reads.append((f"read_id_{i+1}", read))
        
    return reads


# Write simulated reads to FASTQ format
def write_reads_to_fastq(reads, filename):
    
    # Write each read in FASTQ format (4 lines per read)
    with open(filename, 'w') as f:
        for read_id, seq in reads:
            f.write(f"@{read_id}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write('I' * 100 + '\n')        

# Call variants using minimap2 and bcftools pipeline (fully piped)
def call_variants(ref, reads, out_prefix):
    
    # Define output VCF filename
    vcf = f"{out_prefix}.vcf"    
    # Align reads to reference using minimap2
    minimap_proc = subprocess.Popen(["minimap2", "-ax", "sr", ref, reads], stdout=subprocess.PIPE, text=True)    
    # Sort alignments with samtools
    samtools_sort_proc = subprocess.Popen(["samtools", "sort", "-o", "-", "-"], stdin=minimap_proc.stdout, stdout=subprocess.PIPE)
    minimap_proc.stdout.close()
    # Generate pileup with bcftools mpileup
    mpileup_proc = subprocess.Popen(["bcftools", "mpileup", "-f", ref, "-"], stdin=samtools_sort_proc.stdout, stdout=subprocess.PIPE)
    samtools_sort_proc.stdout.close()
    
    # Call variants with bcftools call and write to VCF file
    with open(vcf, "w") as vcf_out:
        subprocess.run(["bcftools", "call", "-mv", "-Ov"], stdin=mpileup_proc.stdout, stdout=vcf_out, check=True)
    mpileup_proc.stdout.close()
    
    return vcf


# Parse VCF file and extract variant information
def parse_vcf(vcf_file):
    
    # List to store parsed variants
    variants = []
    
    # Read VCF file line by line
    with open(vcf_file) as f:
        for line in f:
            # Skip header lines
            if line.startswith("#"):
                continue
            # Extract position, reference, and alternate alleles
            cols = line.strip().split("\t")
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            variants.append((pos, ref, alt))
            
    return variants


# Parse VCF into dictionary for efficient position lookup
def parse_vcf_to_dict(vcf_file):
    
    # Dictionary with position as key, list of (ref, alt) tuples as value
    vcf_dict = defaultdict(list)
    called_variants = parse_vcf(vcf_file)
    
    # Group variants by position
    for pos, ref, alt in called_variants:
        vcf_dict[pos].append((ref, alt))
        
    return vcf_dict


# Compare true mutations against VCF calls and identify matches
def compare_variants(genome_string, indels, vcf_file):
    
    # Load VCF calls into dictionary for efficient lookup
    vcf_lookup = parse_vcf_to_dict(vcf_file)
    comparison_records = []
    true_positives = 0
    
    # Check each true mutation against VCF calls
    for type, pos_0b, identifier, change in indels:
        acc_ref, acc_alt, var_type, pos_1b_vcf_anchor = "", "", "", 0

        # Format SNPs (no anchor base needed)
        if type == "SNP":
            acc_ref = identifier
            acc_alt = change
            var_type = "SNP"
            pos_1b_vcf_anchor = pos_0b + 1
        
        # Format insertions (VCF requires anchor base before insertion)    
        elif type == "INS":
            vcf_anchor_pos_0b = pos_0b - 1
            if vcf_anchor_pos_0b < 0:
                vcf_anchor_pos_0b = pos_0b 
            anchor_base = genome_string[vcf_anchor_pos_0b]
            pos_1b_vcf_anchor = vcf_anchor_pos_0b + 1
            acc_ref = anchor_base
            acc_alt = anchor_base + change
            var_type = "INS"
        
        # Format deletions (VCF requires anchor base before deletion)    
        elif type == "DEL":
            vcf_anchor_pos_0b = pos_0b - 1
            if vcf_anchor_pos_0b < 0:
                vcf_anchor_pos_0b = pos_0b
            anchor_base = genome_string[vcf_anchor_pos_0b]
            pos_1b_vcf_anchor = vcf_anchor_pos_0b + 1
            acc_ref = anchor_base + change
            acc_alt = anchor_base
            var_type = "DEL"
        else:
            continue
        
        # Initialize VCF fields as not found    
        vcf_ref = "---"
        vcf_alt = "---"
        status = "NO"
        vcf_alleles_found = "---"

        # Check if VCF has any calls at this position
        called_at_pos = vcf_lookup.get(pos_1b_vcf_anchor, [])
        
        # If VCF called something at this position
        if called_at_pos:
            vcf_alleles_found = ";".join([f"{r}->{a}" for r, a in called_at_pos])

            # Check if VCF call exactly matches our mutation (True Positive)
            if (acc_ref, acc_alt) in called_at_pos:
                vcf_ref = acc_ref
                vcf_alt = acc_alt
                status = "YES" # True Positive (TP)
                true_positives += 1
        
        # Record comparison result        
        comparison_records.append({
            "ACC_POS(1b)": pos_1b_vcf_anchor,
            "VCF_POS(1b)": pos_1b_vcf_anchor,
            "ACC_REF": acc_ref,
            "VCF_REF": vcf_ref,
            "ACC_ALT": acc_alt,
            "VCF_ALT": vcf_alt,
            "TYPE": var_type,
            "MATCH": status,
            "VCF_ALLELES_FOUND": vcf_alleles_found
        })
        
    return comparison_records, true_positives, vcf_lookup


# Calculate precision, recall, and other performance metrics
def calculate_metrics(indels, true_positives, vcf_lookup): 
    
    # Count total mutations and variants
    total_ground_truth = len(indels)
    total_called_variants = sum(len(v) for v in vcf_lookup.values())
    false_negatives = total_ground_truth - true_positives
    false_positives = total_called_variants - true_positives    
    # Calculate precision (what fraction of calls are correct)
    precision = true_positives / total_called_variants if total_called_variants > 0 else 0
    # Calculate recall (what fraction of true variants were found)
    recall = true_positives / total_ground_truth if total_ground_truth > 0 else 0
    
    return {
        "Total_True_Mutations": total_ground_truth,
        "Total_Called_Variants": total_called_variants,
        "True_Positives": true_positives,
        "False_Negatives": false_negatives,
        "False_Positives": false_positives,
        "Precision": precision,
        "Recall": recall}


# Write comparison results to CSV file
def write_comparison_to_csv(comparison_records, filename):
    
    # Check if there are records to write
    if not comparison_records:
        return

    # Get column names from first record
    fieldnames = comparison_records[0].keys()
    
    # Write all records to CSV
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(comparison_records)

        
# Main function to score variant calls against truth set
def score_variants(genome_string, indels, vcf_file, filename_prefix):
    
    # Compare true mutations with VCF calls
    comparison_records, true_positives, vcf_lookup = compare_variants(genome_string, indels, vcf_file)
    # Save comparison to CSV
    comparison_filename = f"{filename_prefix}_comparison.csv"
    write_comparison_to_csv(comparison_records, comparison_filename)
    # Calculate performance metrics
    metrics = calculate_metrics(indels, true_positives, vcf_lookup)
    
    return metrics