#!/bin/bash

# Paths
#BUSCO_IDS_FILE="not192busco_counts.tsv"  # File containing BUSCO IDs (first column)

BUSCO_IDS_FILE="all_buscos.txt"  # File containing BUSCO IDs (first column)

BUSCO_OUTPUT_FOLDER="/mnt/data/dtol_review/busco_out_AllGenomes"  # BUSCO directories
#OUTPUT_DIR="grouped_busco_fastas/secondround_withallgenes"  # Directory for grouped FASTA files
OUTPUT_DIR="grouped_busco_fastas/review_allbuscos_442species"

# Create the output directory
mkdir -p "$OUTPUT_DIR"

# Read BUSCO IDs from the file
while IFS=$'\t' read -r busco_id _; do
    # Skip the header line
    [[ "$busco_id" == "gene_id" ]] && continue

    echo "Processing BUSCO ID: $busco_id"

    # Output file for this BUSCO ID
    output_fasta="$OUTPUT_DIR/$busco_id.faa"
    # Check if the output file already exists
    if [[ -f "$output_fasta" ]]; then
    	echo "File $output_fasta already exists. Skipping..."
	continue
    fi

    # Clear or create the output FASTA file for the current BUSCO ID
    > "$output_fasta"

    # Find all matching .faa files for this BUSCO ID
    find "$BUSCO_OUTPUT_FOLDER" -type f -path "*/run_eukaryota_odb10/busco_sequences/*" -name "$busco_id.faa" | while read -r fasta_file; do
        # Extract the species identifier (third-to-last directory level)
        species_id=$(basename "$(dirname "$(dirname "$(dirname "$(dirname "$fasta_file")")")")")

        # Print the species ID and busco ID
        echo "Adding sequence from species: $species_id for BUSCO ID: $busco_id"

        # Use Biopython to process the first sequence and modify the header
        python3 - <<EOF
from Bio import SeqIO

# Open the input FASTA file
with open("$fasta_file", "r") as handle:
    # Parse the FASTA file and get the first sequence
    for i, record in enumerate(SeqIO.parse(handle, "fasta")):
        if i == 0:  # Process only the first sequence
            # Modify the header by prepending species ID and the original FASTA header
            record.id = "$species_id|" + record.id
            record.description = ""  # Remove any description
            
            # Write the modified sequence to the output FASTA file
            with open("$output_fasta", "a") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
            break  # Stop after the first sequence is processed
EOF
    done

    echo "Finished BUSCO ID: $busco_id, saved to $output_fasta"
done < "$BUSCO_IDS_FILE"

echo "All grouped sequences saved in $OUTPUT_DIR"

