import os
import pandas as pd

def find_absent_buscos(root_dir, top_buscos):
    """
    Find the species where the most highly occurring BUSCO genes are absent.
    
    Args:
        root_dir (str): Path to the root directory containing BUSCO result subdirectories.
        top_buscos (list): List of BUSCO IDs to check for absence.

    Returns:
        pd.DataFrame: A DataFrame showing species absence for each BUSCO ID.
    """
    busco_absence = {busco_id: [] for busco_id in top_buscos}
    total_dirs = len(os.listdir(root_dir))
    processed_dirs = 0

    # Process each species directory
    for species_dir in os.listdir(root_dir):
        species_path = os.path.join(root_dir, species_dir)
        if os.path.isdir(species_path):  # Check if it's a directory
            processed_dirs += 1
            print(f"Processing directory {processed_dirs}/{total_dirs}: {species_dir}")

            # Path to full_table.tsv
            busco_run_path = os.path.join(species_path, "run_eukaryota_odb10")
            full_table_path = os.path.join(busco_run_path, "full_table.tsv")

            if os.path.exists(full_table_path):
                try:
                    # Read the file, skipping lines starting with '#'
                    data = pd.read_csv(full_table_path, sep="\t", comment="#", header=None)

                    # Keep BUSCOs marked as "Complete" OR "Duplicated"
                    present_buscos = set(
                        data[data[1].isin(["Complete", "Duplicated"])][0].tolist()
                    )

                    # Check for absence of top BUSCO IDs
                    for busco_id in top_buscos:
                        if busco_id not in present_buscos:
                            busco_absence[busco_id].append(species_dir)
                except Exception as e:
                    print(f"Error reading {full_table_path}: {e}")
            else:
                print(f"Warning: 'full_table.tsv' not found in {busco_run_path}")

    # Create a DataFrame from the absence dictionary
    busco_absence_df = pd.DataFrame.from_dict(busco_absence, orient="index").transpose()
    busco_absence_df.columns = top_buscos
    busco_absence_df.fillna("", inplace=True)  # Replace NaNs with empty strings
    return busco_absence_df

# Example usage
root_dir = "/mnt/data/dtol_review/busco_out_AllGenomes"
top_buscos = ["1314980at2759", "1324510at2759", "290630at2759", "299766at2759",
              "320059at2759", "901894at2759", "905026at2759", "937275at2759",
              "166920at2759", "674160at2759", "324863at2759", "330169at2759"]

# Get absences for top BUSCOs
busco_absence_df = find_absent_buscos(root_dir, top_buscos)

# Save the results to a TSV file
output_file = os.path.join(root_dir, "busco_absences_duplicated_notincluded.tsv")
busco_absence_df.to_csv(output_file, sep="\t", index=False)
print(f"BUSCO absences saved to: {output_file}")

