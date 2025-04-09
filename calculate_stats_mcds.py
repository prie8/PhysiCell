import os
import glob
import numpy as np
import pandas as pd
from pyMCDS import pyMCDS  # Requires pyMCDS, numpy, pandas

# --- Configuration ---
# Set the path to your PhysiCell output directory
output_directory = (
    "./output/"  # Should be relative or absolute path to the output folder
)

# Cell type IDs used in the simulation (Match your cell definitions in XML)
biofilm_type_id = 0
probiotic_type_id = 1

# Custom data variable name for IL10 secretion state (check pyMCDS output/cell_df.columns if needed)
il10_custom_data_var = "custom:secreting_il10"
# ---------------------

print(f"Analyzing output files in: {output_directory}")
print(f"Determining initial counts from first file, final state counts from last file.")

# --- Find simulation output files ---
file_pattern = os.path.join(output_directory, "output*.xml")
# Use sorted list of FULL PATHS first
xml_files_full_path = sorted(glob.glob(file_pattern))

if not xml_files_full_path:
    print(f"Error: No output*.xml files found in '{output_directory}'.")
    print("Please ensure the output directory is correct and the simulation ran.")
    exit()

# --- Process FIRST file to get initial counts ---
initial_biofilm_count = 0.0  # Use float for division later
initial_probiotic_count = 0.0
first_xml_full_path = xml_files_full_path[0]
first_xml_basename = os.path.basename(first_xml_full_path)  # Extract just the filename
print(f"\nReading initial counts from: {first_xml_basename} in {output_directory}")

try:
    # Load the first simulation state using basename and directory path
    mcds_init = pyMCDS(first_xml_basename, output_directory)  # CORRECTED CALL
    # Get the initial cell data frame
    cell_df_init = mcds_init.get_cell_df()

    if not cell_df_init.empty:
        # Count cells of each type
        counts_init = cell_df_init["cell_type"].value_counts().to_dict()
        initial_biofilm_count = float(counts_init.get(float(biofilm_type_id), 0))
        initial_probiotic_count = float(counts_init.get(float(probiotic_type_id), 0))
    else:
        print(f"Warning: No cells found in the initial file {first_xml_basename}.")
except Exception as e:
    print(f"\nError processing initial file {first_xml_basename}: {e}")
    exit()


# --- Process LAST file to get final counts and states ---
final_biofilm_live = 0
final_biofilm_dead = 0
final_probiotic_magenta = 0
final_probiotic_yellow = 0
final_probiotic_cyan = 0
final_probiotic_green = 0
final_probiotic_dead = 0
final_total_cells = 0
last_xml_full_path = xml_files_full_path[-1]
last_xml_basename = os.path.basename(last_xml_full_path)  # Extract just the filename
print(
    f"\nReading final counts/states from:   {last_xml_basename} in {output_directory}"
)

try:
    # Load the last simulation state using basename and directory path
    mcds_final = pyMCDS(last_xml_basename, output_directory)  # CORRECTED CALL
    # Get the final cell data frame
    cell_df_final = mcds_final.get_cell_df()

    if not cell_df_final.empty:
        final_total_cells = len(cell_df_final)
        biofilm_cells = cell_df_final[
            cell_df_final["cell_type"] == float(biofilm_type_id)
        ]
        probiotic_cells = cell_df_final[
            cell_df_final["cell_type"] == float(probiotic_type_id)
        ]

        # Count biofilm states
        final_biofilm_live = biofilm_cells[biofilm_cells["dead"] == 0].shape[0]
        final_biofilm_dead = biofilm_cells[biofilm_cells["dead"] == 1].shape[0]

        # Count probiotic states (logic unchanged)
        final_probiotic_dead = probiotic_cells[probiotic_cells["dead"] == 1].shape[0]
        live_probiotics = probiotic_cells[probiotic_cells["dead"] == 0].copy()
        if not live_probiotics.empty:
            if il10_custom_data_var in live_probiotics.columns:
                live_probiotics["secreting_dnase"] = True
                live_probiotics["secreting_il10"] = (
                    live_probiotics[il10_custom_data_var] > 0.5
                )
                final_probiotic_magenta = live_probiotics[
                    live_probiotics["secreting_dnase"]
                    & live_probiotics["secreting_il10"]
                ].shape[0]
                final_probiotic_yellow = live_probiotics[
                    live_probiotics["secreting_dnase"]
                    & ~live_probiotics["secreting_il10"]
                ].shape[0]
                final_probiotic_cyan = live_probiotics[
                    ~live_probiotics["secreting_dnase"]
                    & live_probiotics["secreting_il10"]
                ].shape[0]
                final_probiotic_green = live_probiotics[
                    ~live_probiotics["secreting_dnase"]
                    & ~live_probiotics["secreting_il10"]
                ].shape[0]
            else:
                print(
                    f"Warning: Custom data column '{il10_custom_data_var}' not found."
                )
                final_probiotic_green = live_probiotics.shape[0]
    else:
        print(f"Warning: No cells found in the final file {last_xml_basename}.")

except Exception as e:
    print(f"\nError processing final file {last_xml_basename}: {e}")
    exit()


# --- Calculate Percentages based on FINAL LIVE vs DEAD biofilm ---
final_biofilm_total = final_biofilm_live + final_biofilm_dead
if final_biofilm_total > 0:
    percent_biofilm_live_final = (final_biofilm_live / final_biofilm_total) * 100.0
    percent_biofilm_dead_final = (final_biofilm_dead / final_biofilm_total) * 100.0
else:
    percent_biofilm_live_final = 0.0 if final_biofilm_live == 0 else np.nan
    percent_biofilm_dead_final = 0.0 if final_biofilm_dead == 0 else np.nan


# --- Print Results ---
print("\n--- Simulation Final Statistics ---")
# ... (rest of print statements unchanged) ...
print(f"Initial Biofilm Cells:      {int(initial_biofilm_count)}")
print(f"Initial Probiotic Cells:    {int(initial_probiotic_count)}")
print(
    f"Initial Total Cells:        {int(initial_biofilm_count + initial_probiotic_count)}"
)
print("-" * 30)
print(f"Final Biofilm Cells (Live): {final_biofilm_live}")
print(f"Final Biofilm Cells (Dead): {final_biofilm_dead}")
print(f"  (Total Final Biofilm:    {final_biofilm_total})")
print("-" * 30)
print(f"Final Probiotic Cells (Magenta): {final_probiotic_magenta}")
print(f"Final Probiotic Cells (Yellow):  {final_probiotic_yellow}")
print(f"Final Probiotic Cells (Cyan):    {final_probiotic_cyan}")
print(f"Final Probiotic Cells (Green):   {final_probiotic_green}")
print(f"Final Probiotic Cells (Dead):  {final_probiotic_dead}")
print(
    f"  (Sum Probiotic Live: {final_probiotic_magenta + final_probiotic_yellow + final_probiotic_cyan + final_probiotic_green})"
)
print("-" * 30)
print(f"Final Total Cells:          {final_total_cells}")
print("-" * 30)
if not np.isnan(percent_biofilm_dead_final):
    print(f"Biofilm Live Pct (Live/Total Final): {percent_biofilm_live_final:.2f}%")
    print(f"Biofilm Reduction (Dead/Total Final):{percent_biofilm_dead_final:.2f}%")
else:
    print("Biofilm percentages undefined (no biofilm cells found in final state).")
print("-" * 30)

print("Done.")
