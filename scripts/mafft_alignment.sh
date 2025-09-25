#!/bin/bash

# Input and output directories
INPUT_DIR="grouped_busco_fastas/review_allbuscos_442species/"  # Current directory with FASTA files
OUTPUT_DIR="grouped_busco_fastas/review_allbuscos_442species/msas"  # Directory to save aligned files

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .faa files in the input directory
for fasta_file in "$INPUT_DIR"/*.faa; do
    # Extract the base name without the extension
    base_name=$(basename "$fasta_file" .faa)

    # Define the output file path
    output_file="$OUTPUT_DIR/$base_name.msa"

    echo "Aligning $fasta_file to $output_file"

    # Run MAFFT with specified options
    mafft --auto --thread 10 --maxiterate 1000 "$fasta_file" > "$output_file"

    # Check if MAFFT succeeded
    if [[ $? -eq 0 ]]; then
        echo "Successfully aligned $fasta_file"
    else
        echo "Error aligning $fasta_file" >&2
    fi
done

echo "All alignments completed. Results are in $OUTPUT_DIR"

