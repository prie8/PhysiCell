import os
import glob
import numpy as np
import pandas as pd
from pyMCDS import pyMCDS # Requires pyMCDS, numpy, pandas

# --- Configuration ---
# Set the path to your PhysiCell output directory
output_directory = './output/' # Should be relative path from where script is run, or absolute path

# Cell type IDs used in the simulation (Match your cell definitions in XML)
biofilm_type_id = 0
probiotic_type_id = 1

# Custom data variable name for IL10 secretion state (check pyMCDS output/cell_df.columns if needed)
il10_custom_data_var = 'custom:secreting_il10'

# Output CSV filename
output_csv_file = os.path.join(output_directory,'simulation_stats.csv')
# ---------------------

print(f"Analyzing output files in: {output_directory}")
print(f"Determining initial counts from first file, final state counts from last file.")

# --- Find simulation output files ---
file_pattern = os.path.join(output_directory, 'output*.xml')
xml_files = sorted(glob.glob(file_pattern))

if not xml_files:
    print(f"Error: No output*.xml files found in '{output_directory}'.")
    print("Please ensure the output directory is correct and the simulation ran.")
    exit()

print(f"Found {len(xml_files)} simulation time steps to process...")

# --- Process FIRST file to get initial counts ---
initial_biofilm_count = 0.0 # Use float for division later
initial_probiotic_count = 0.0
first_xml_file = xml_files[0]
print(f"\nReading initial counts from: {os.path.basename(first_xml_file)}")

try:
    mcds_init = pyMCDS(os.path.basename(first_xml_file), output_directory) # Use basename
    cell_df_init = mcds_init.get_cell_df()
    if not cell_df_init.empty:
        counts_init = cell_df_init['cell_type'].value_counts().to_dict()
        initial_biofilm_count = float(counts_init.get(float(biofilm_type_id), 0))
        initial_probiotic_count = float(counts_init.get(float(probiotic_type_id), 0))
    else:
         print(f"Warning: No cells found in the initial file {first_xml_file}.")
except Exception as e:
    print(f"\nError processing initial file {first_xml_file}: {e}")
    exit()


# --- Process All Files for CSV Output ---
all_stats_data = [] # List to store stats dictionaries from each time step

for i, xml_full_path in enumerate(xml_files):
    xml_basename = os.path.basename(xml_full_path)
    if (i == 0) or ((i + 1) % 10 == 0) or (i + 1 == len(xml_files)): # Print progress
        print(f"  Processing file {i+1}/{len(xml_files)}: {xml_basename}")

    try:
        mcds = pyMCDS(xml_basename, output_directory) # Use basename
        current_time = mcds.get_time()
        cell_df = mcds.get_cell_df()

        # Initialize counts for this time step
        biofilm_live = 0; biofilm_dead = 0
        probiotic_magenta = 0; probiotic_yellow = 0; probiotic_cyan = 0; probiotic_green = 0; probiotic_dead = 0
        total_cells = 0

        if not cell_df.empty:
            total_cells = len(cell_df)
            biofilm_cells = cell_df[cell_df['cell_type'] == float(biofilm_type_id)]
            probiotic_cells = cell_df[cell_df['cell_type'] == float(probiotic_type_id)]

            # Count biofilm states
            biofilm_live = biofilm_cells[biofilm_cells['dead'] == 0].shape[0]
            biofilm_dead = biofilm_cells[biofilm_cells['dead'] == 1].shape[0]

            # Count probiotic states
            probiotic_dead = probiotic_cells[probiotic_cells['dead'] == 1].shape[0]
            live_probiotics = probiotic_cells[probiotic_cells['dead'] == 0].copy()
            if not live_probiotics.empty:
                if il10_custom_data_var in live_probiotics.columns:
                    live_probiotics['secreting_dnase'] = True # Assumed TRUE
                    live_probiotics['secreting_il10'] = live_probiotics[il10_custom_data_var] > 0.5
                    probiotic_magenta = live_probiotics[live_probiotics['secreting_dnase'] & live_probiotics['secreting_il10']].shape[0]
                    probiotic_yellow = live_probiotics[live_probiotics['secreting_dnase'] & ~live_probiotics['secreting_il10']].shape[0]
                    probiotic_cyan = live_probiotics[~live_probiotics['secreting_dnase'] & live_probiotics['secreting_il10']].shape[0]
                    probiotic_green = live_probiotics[~live_probiotics['secreting_dnase'] & ~live_probiotics['secreting_il10']].shape[0]
                else:
                    probiotic_yellow = live_probiotics.shape[0] # Default if custom data missing
                    # if i == 0: print(f"Warning: Custom data column '{il10_custom_data_var}' not found.")
        else:
             if i == 0: print(f"Warning: No cells found in file {xml_basename}.")

        # Store stats for CSV
        stats_dict = {
            'time_min': current_time, 'total_cells': total_cells,
            'biofilm_live': biofilm_live, 'biofilm_dead': biofilm_dead,
            'probiotic_live_magenta': probiotic_magenta, 'probiotic_live_yellow': probiotic_yellow,
            'probiotic_live_cyan': probiotic_cyan, 'probiotic_live_green': probiotic_green,
            'probiotic_dead': probiotic_dead
        }
        all_stats_data.append(stats_dict)

        # --- Store final counts explicitly for summary printout ---
        if i == len(xml_files) - 1:
            final_biofilm_live = biofilm_live
            final_biofilm_dead = biofilm_dead
            final_probiotic_magenta = probiotic_magenta
            final_probiotic_yellow = probiotic_yellow
            final_probiotic_cyan = probiotic_cyan
            final_probiotic_green = probiotic_green
            final_probiotic_dead = probiotic_dead
            final_total_cells = total_cells
        # --- End storing final counts ---

    except Exception as e:
        print(f"\nError processing file {xml_basename}: {e}")
        # Continue processing other files if possible

# --- Create DataFrame and Save to CSV ---
if all_stats_data:
    stats_df = pd.DataFrame(all_stats_data)
    # Add percentage columns to DataFrame
    if initial_biofilm_count > 0:
        stats_df['biofilm_pct_reduction_vs_initial'] = 100.0 * (initial_biofilm_count - stats_df['biofilm_live']) / initial_biofilm_count
    else: stats_df['biofilm_pct_reduction_vs_initial'] = np.nan
    stats_df['biofilm_total_final_step'] = stats_df['biofilm_live'] + stats_df['biofilm_dead']
    stats_df['biofilm_pct_reduction_vs_final_step'] = np.where( stats_df['biofilm_total_final_step'] > 0, 100.0 * stats_df['biofilm_dead'] / stats_df['biofilm_total_final_step'], 0.0)
    # Define column order
    column_order = [
        'time_min', 'total_cells', 'biofilm_live', 'biofilm_dead', 'biofilm_total_final_step',
        'biofilm_pct_reduction_vs_initial', 'biofilm_pct_reduction_vs_final_step',
        'probiotic_live_magenta', 'probiotic_live_yellow', 'probiotic_live_cyan',
        'probiotic_live_green', 'probiotic_dead' ]
    stats_df = stats_df.reindex(columns=column_order, fill_value=0)
    # Save CSV
    try:
        stats_df.to_csv(output_csv_file, index=False, float_format='%.4f')
        print(f"\nSuccessfully saved time series stats to: {output_csv_file}")
    except Exception as e:
        print(f"\nError saving CSV file: {e}")
else:
    print("\nNo data processed, CSV file not saved.")


# --- Calculate FINAL Percentages using final counts ---
final_biofilm_total = final_biofilm_live + final_biofilm_dead
if final_biofilm_total > 0:
    percent_biofilm_live_final = (final_biofilm_live / final_biofilm_total) * 100.0
    percent_biofilm_dead_final = (final_biofilm_dead / final_biofilm_total) * 100.0
else:
    percent_biofilm_live_final = 0.0 if final_biofilm_live == 0 else np.nan
    percent_biofilm_dead_final = 0.0 if final_biofilm_dead == 0 else np.nan

# --- ADDED BACK: Print Final Summary to Console ---
print("\n--- Simulation Final Statistics ---")
print(f"Initial Biofilm Cells:      {int(initial_biofilm_count)}")
print(f"Initial Probiotic Cells:    {int(initial_probiotic_count)}")
print(f"Initial Total Cells:        {int(initial_biofilm_count + initial_probiotic_count)}")
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
probiotic_live_sum = final_probiotic_magenta + final_probiotic_yellow + final_probiotic_cyan + final_probiotic_green
print(f"  (Sum Probiotic Live: {probiotic_live_sum})")
print("-" * 30)
print(f"Final Total Cells:          {final_total_cells}")
print("-" * 30)
if not np.isnan(percent_biofilm_dead_final):
    print(f"Biofilm Live Pct (Live/Total Final): {percent_biofilm_live_final:.2f}%")
    print(f"Biofilm Reduction (Dead/Total Final):{percent_biofilm_dead_final:.2f}%")
else:
     print("Biofilm percentages undefined (no biofilm cells found in final state).")
# Optionally calculate and print reduction vs initial
if initial_biofilm_count > 0:
     percent_biofilm_killed_vs_initial = 100.0 * (initial_biofilm_count - final_biofilm_live) / initial_biofilm_count
     print(f"Biofilm Reduction (vs Initial):    {percent_biofilm_killed_vs_initial:.2f}%")

print("-" * 30)
print("Done.")
