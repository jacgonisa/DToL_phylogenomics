import os
import pandas as pd

def count_complete_buscos(busco_output_folder):
    """
    Count the occurrences of "Complete" or "Duplicated" BUSCOs across all species directories.
    Ensure each BUSCO is counted only once per species.

    Args:
        busco_output_folder (str): Path to the BUSCO output folder containing species directories.

    Returns:
        dict: Dictionary with BUSCO IDs as keys and counts as values.
    """
    busco_counts = {}

    # Iterate through directories in the BUSCO output folder
    for root, dirs, files in os.walk(busco_output_folder):
        for directory in dirs:
            print(f"Processing directory: {directory}")
            full_table_path = os.path.join(root, directory, "run_eukaryota_odb10", "full_table.tsv")

            # Check if the full_table.tsv file exists
            if os.path.exists(full_table_path):
                # Read the file, ignoring lines starting with '#'
                data = pd.read_csv(
                    full_table_path,
                    sep="\t",
                    comment="#",
                    header=None,
                    usecols=[0, 1],  # Only need columns 0 (BUSCO ID) and 1 (status)
                )
                
                # Filter rows for "Complete" or "Duplicated"
                valid_buscos = data[data[1].isin(["Complete", "Duplicated"])][0].drop_duplicates()
                
                # Update the counts for each BUSCO ID
                for busco in valid_buscos:
                    if busco in busco_counts:
                        busco_counts[busco] += 1
                    else:
                        busco_counts[busco] = 1

    return busco_counts

# Example usage
busco_output_folder = "/mnt/data/dtol_review/busco_out_AllGenomes"
busco_counts = count_complete_buscos(busco_output_folder)

# Save results to a TSV file
output_file = "busco_counts.tsv"
with open(output_file, "w") as f:
    f.write("BUSCO_ID\tOccurrences\n")
    for busco_id, count in sorted(busco_counts.items(), key=lambda x: -x[1]):  # Sort by count (descending)
        f.write(f"{busco_id}\t{count}\n")

print(f"BUSCO counts saved to {output_file}")

