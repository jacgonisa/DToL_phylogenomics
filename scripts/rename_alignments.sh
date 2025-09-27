#!/bin/bash

# Directory containing the .msa files
MSA_DIR="grouped_busco_fastas/review_allbuscos_442species/trimmed_msas/"

# Loop through all .msa files in the directory
for msa_file in "$MSA_DIR"/*.msa; do
	    # Generate the new file name
	        renamed_file="${msa_file%.msa}.renamed.msa"
		    
		    # Process the file
		        awk '/^>/ {split($0, arr, "\\|"); print arr[1]} !/^>/ {print}' "$msa_file" > "$renamed_file"
		done

		echo "Renamed sequence headers saved to new .msa files with .renamed.msa extension."

