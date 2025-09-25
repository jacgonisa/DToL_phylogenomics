#!/bin/bash

# Input and output directories
INPUT_DIR="grouped_busco_fastas/review_allbuscos_442species/msas/"  # Directory containing .msa files
OUTPUT_DIR="grouped_busco_fastas/review_allbuscos_442species/trimmed_msas/"  # Directory to save trimmed alignments

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .msa files in the input directory
for msa_file in "$INPUT_DIR"/*.msa; do
    # Extract the base name without the extension
    base_name=$(basename "$msa_file" .msa)

    # Define the output file path
    output_file="$OUTPUT_DIR/$base_name.trimmed.msa"

    echo "Trimming $msa_file to $output_file"

    # Run trimAl with the "gappyout" option
    trimal -in "$msa_file" -out "$output_file" -gappyout

    # Check if trimAl succeeded
    if [[ $? -eq 0 ]]; then
        echo "Successfully trimmed $msa_file"
    else
        echo "Error trimming $msa_file" >&2
    fi
done

echo "All trimming completed. Results are in $OUTPUT_DIR"


