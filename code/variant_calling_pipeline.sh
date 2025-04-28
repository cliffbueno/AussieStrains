#!/bin/bash

# Set paths to tools and files
REFERENCE="/scratch/alpine/clbd1748/Australia/BradyStrains/BradyReference/BradyReference.fasta"  # Ref genome
GENOME_DIR="/scratch/alpine/clbd1748/Australia/BradyStrains/"  # 53 strains
OUTPUT_DIR="/scratch/alpine/clbd1748/Australia/variant_calling/"  # Output directory
LOG_FILE="variant_calling.log"      # Log file for output text
THREADS=32                   # Number of threads to use

# Step 0: Create output directory and log file
mkdir -p $OUTPUT_DIR
echo "Starting variant calling pipeline..." > $LOG_FILE

# Step 1: Index the reference genome
echo "Indexing reference genome..." | tee -a $LOG_FILE
bwa index $REFERENCE 2>&1 | tee -a $LOG_FILE
samtools faidx $REFERENCE 2>&1 | tee -a $LOG_FILE
gatk CreateSequenceDictionary -R $REFERENCE -O "${REFERENCE%.fasta}.dict" 2>&1 | tee -a $LOG_FILE

# Step 2: Loop through genomes and align them
echo "Aligning genomes to the reference..." | tee -a $LOG_FILE
for GENOME in $GENOME_DIR*.fasta; do
    BASENAME=$(basename $GENOME .fasta)  # Extract file name without extension
    echo "Processing $BASENAME..." | tee -a $LOG_FILE

    # Align genome to reference
    bwa mem -t $THREADS $REFERENCE $GENOME > $OUTPUT_DIR${BASENAME}.sam 2>> $LOG_FILE

    # Convert SAM to BAM, sort, and index
    samtools view -@ $THREADS -S -b $OUTPUT_DIR${BASENAME}.sam | \
        samtools sort -@ $THREADS -o $OUTPUT_DIR${BASENAME}_sorted.bam 2>> $LOG_FILE
    samtools index $OUTPUT_DIR${BASENAME}_sorted.bam 2>> $LOG_FILE

    # Clean up SAM file
    rm $OUTPUT_DIR${BASENAME}.sam
done

# Step 3: Call variants using bcftools
echo "Calling variants with bcftools..." | tee -a $LOG_FILE
BAM_FILES=$(ls $OUTPUT_DIR*_sorted.bam | tr '\n' ' ')  # List all sorted BAM files
bcftools mpileup -Ou -f $REFERENCE $BAM_FILES 2>> $LOG_FILE | \
    bcftools call -mv -Ov -o $OUTPUT_DIR/variants_raw.vcf 2>> $LOG_FILE

# Step 4: Filter the variants
echo "Filtering variants..." | tee -a $LOG_FILE
bcftools filter -i 'QUAL > 30 && DP > 10' $OUTPUT_DIR/variants_raw.vcf -o $OUTPUT_DIR/variants_filtered.vcf 2>> $LOG_FILE

# Final output
echo "Variant calling pipeline complete. Results are in the directory: $OUTPUT_DIR" | tee -a $LOG_FILE
