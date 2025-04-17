import os
import glob
import numpy as np
import pandas as pd
from pyMCDS import pyMCDS
import sys
import matplotlib.pyplot as plt # Import matplotlib

# --- Configuration ---
# Set the path to your PhysiCell output directory (Updated based on XML)
output_directory = './output_P1/'
if len(sys.argv) > 1: output_directory = sys.argv[1]

# Cell type IDs
biofilm_type_id = 0
probiotic_type_id = 1

# --- Custom data variable name for IL10 state (No longer primary method) ---
# il10_custom_data_var = 'custom:secreting_il10' # Kept for reference / potential fallback

# Output files
output_csv_file = 'simulation_stats.csv'
output_plot_file = 'simulation_summary_plots.png'

# Target reduction range for plotting
target_reduction_min = 50.0
target_reduction_max = 70.0
# ---------------------

print(f"Analyzing output files in: {output_directory}")
print(f"Calculating stats and generating plots...")

# --- Find simulation output files ---
file_pattern = os.path.join(output_directory, 'output*.xml')
xml_files = sorted(glob.glob(file_pattern))

if not xml_files:
    print(f"Error: No output*.xml files found in '{output_directory}'.")
    exit()

print(f"Found {len(xml_files)} simulation time steps to process...")

# --- Process All Files ---
all_stats_data = []
initial_biofilm_count = 0.0
initial_probiotic_count = 0.0

for i, xml_full_path in enumerate(xml_files):
    xml_basename = os.path.basename(xml_full_path)
    # Print progress less often for faster processing
    if (i == 0) or ((i + 1) % 50 == 0) or (i + 1 == len(xml_files)):
        print(f"  Processing file {i+1}/{len(xml_files)}: {xml_basename}")

    try:
        mcds = pyMCDS(xml_basename, output_directory)
        current_time = mcds.get_time()
        cell_df = mcds.get_cell_df()

        # ADD THIS: Print columns only for the first file for verification
        # if i == 0 and not cell_df.empty:
        #     print("\nDataFrame Columns Found:")
        #     print(list(cell_df.columns))
        #     print("-" * 30)
        # END ADD

        # Initialize counts
        biofilm_live = 0; biofilm_dead = 0
        probiotic_lime = 0; probiotic_yellow = 0; probiotic_cyan = 0; probiotic_magenta = 0; probiotic_dead = 0
        total_cells = 0

        if not cell_df.empty:
            total_cells = len(cell_df)
            biofilm_cells = cell_df[cell_df['cell_type'] == float(biofilm_type_id)]
            probiotic_cells = cell_df[cell_df['cell_type'] == float(probiotic_type_id)]

            biofilm_live = biofilm_cells[biofilm_cells['dead'] == 0].shape[0]
            biofilm_dead = biofilm_cells[biofilm_cells['dead'] == 1].shape[0]

            probiotic_dead = probiotic_cells[probiotic_cells['dead'] == 1].shape[0]
            live_probiotics = probiotic_cells[probiotic_cells['dead'] == 0].copy()

            if not live_probiotics.empty:
                # --- Determine state using Correct Secretion Rate Column Names ---
                dnase_rate_col = 'secretion_rates_y' # CORRECTED based on user feedback (DNase index 1)
                il10_rate_col = 'secretion_rates_z'  # CORRECTED based on user feedback (IL10 index 2)
                # Use a small threshold > 0 for checking rates
                rate_threshold = 1e-6

                # Check if the necessary rate columns exist
                if dnase_rate_col in live_probiotics.columns and il10_rate_col in live_probiotics.columns:

                    # Determine boolean states based on rates
                    live_probiotics['is_secreting_dnase'] = live_probiotics[dnase_rate_col] > rate_threshold
                    live_probiotics['is_secreting_il10'] = live_probiotics[il10_rate_col] > rate_threshold

                    # Calculate counts based on C++ color logic using rates
                    # Lime = DNase AND IL10 -> Killing state
                    probiotic_lime = live_probiotics[live_probiotics['is_secreting_dnase'] & live_probiotics['is_secreting_il10']].shape[0]
                    # Yellow = DNase AND NOT IL10 -> Killing state
                    probiotic_yellow = live_probiotics[live_probiotics['is_secreting_dnase'] & ~live_probiotics['is_secreting_il10']].shape[0]
                    # Cyan = NOT DNase AND IL10 -> IL10 Secretion only state
                    probiotic_cyan = live_probiotics[~live_probiotics['is_secreting_dnase'] & live_probiotics['is_secreting_il10']].shape[0]
                    # Magenta = NOT DNase AND NOT IL10 -> Idle state
                    probiotic_magenta = live_probiotics[~live_probiotics['is_secreting_dnase'] & ~live_probiotics['is_secreting_il10']].shape[0]

                else:
                    # Fallback if rate columns are missing
                    probiotic_yellow = live_probiotics.shape[0] # Assign all to a default state
                    probiotic_lime = 0; probiotic_cyan = 0; probiotic_magenta = 0
                    # Updated warning message
                    if i == 0: print(f"Warning: Expected secretion rate columns ('{dnase_rate_col}', '{il10_rate_col}') not found. Using fallback state.")
                # --- END Modification for state calculation ---
        else:
             if i == 0: print(f"Warning: No cells found in file {xml_basename}.")

        # Store initial counts from the very first file (i=0)
        if i == 0:
            initial_biofilm_count = float(biofilm_live + biofilm_dead)
            initial_probiotic_count = float(probiotic_lime + probiotic_yellow + probiotic_cyan + probiotic_magenta + probiotic_dead)

        # Store stats for this time step
        stats_dict = {
            'time_min': current_time, 'total_cells': total_cells,
            'biofilm_live': biofilm_live, 'biofilm_dead': biofilm_dead,
            'probiotic_live_lime': probiotic_lime, 'probiotic_live_yellow': probiotic_yellow,
            'probiotic_live_cyan': probiotic_cyan, 'probiotic_live_magenta': probiotic_magenta,
            'probiotic_dead': probiotic_dead
        }
        all_stats_data.append(stats_dict)

        # Store final counts for summary printout
        if i == len(xml_files) - 1:
            final_biofilm_live = biofilm_live; final_biofilm_dead = biofilm_dead
            final_probiotic_lime = probiotic_lime
            final_probiotic_yellow = probiotic_yellow
            final_probiotic_cyan = probiotic_cyan
            final_probiotic_magenta = probiotic_magenta
            final_probiotic_dead = probiotic_dead
            final_total_cells = total_cells

    except Exception as e:
        print(f"\nError processing file {xml_basename}: {e}")
        # Continue processing other files if possible

# --- Create DataFrame ---
if all_stats_data:
    stats_df = pd.DataFrame(all_stats_data)

    # Calculate derived columns
    stats_df['probiotic_live_total'] = stats_df['probiotic_live_lime'] + stats_df['probiotic_live_yellow'] + \
                                       stats_df['probiotic_live_cyan'] + stats_df['probiotic_live_magenta']
    stats_df['biofilm_total'] = stats_df['biofilm_live'] + stats_df['biofilm_dead']

    if initial_biofilm_count > 0:
        stats_df['biofilm_pct_reduction_vs_initial'] = 100.0 * (initial_biofilm_count - stats_df['biofilm_live']) / initial_biofilm_count
    else:
        stats_df['biofilm_pct_reduction_vs_initial'] = np.nan

    stats_df['biofilm_pct_reduction_vs_final_step'] = np.where(
        stats_df['biofilm_total'] > 0,
        100.0 * stats_df['biofilm_dead'] / stats_df['biofilm_total'],
        0.0
        )

    # Define column order for CSV
    column_order = [
        'time_min', 'total_cells', 'biofilm_live', 'biofilm_dead', 'biofilm_total',
        'biofilm_pct_reduction_vs_initial', 'biofilm_pct_reduction_vs_final_step',
        'probiotic_live_total', 'probiotic_live_lime', 'probiotic_live_yellow',
        'probiotic_live_cyan', 'probiotic_live_magenta', 'probiotic_dead'
    ]
    stats_df = stats_df.reindex(columns=column_order, fill_value=0)

    # --- Save to CSV ---
    try:
        csv_full_path = os.path.join(output_directory, output_csv_file) # Save inside output dir
        stats_df.to_csv(csv_full_path, index=False, float_format='%.4f')
        print(f"\nSuccessfully saved time series stats to: {csv_full_path}")
    except Exception as e:
        print(f"\nError saving CSV file: {e}")

    # --- Generate Plots ---
    print(f"\nGenerating summary plots...")
    try:
        fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True)

        # Plot 1: Cell Counts
        axs[0].plot(stats_df['time_min'], stats_df['biofilm_live'], 'r-', label='Biofilm (Live)')
        axs[0].plot(stats_df['time_min'], stats_df['probiotic_live_total'], 'b-', label='Probiotic (Live)')
        axs[0].plot(stats_df['time_min'], stats_df['total_cells'], 'k--', label='Total Cells')
        axs[0].set_ylabel("Cell Count")
        axs[0].set_title("Cell Populations Over Time")
        axs[0].legend()
        axs[0].grid(True, linestyle=':')

        # Plot 2: Biofilm Reduction
        axs[1].plot(stats_df['time_min'], stats_df['biofilm_pct_reduction_vs_initial'], 'b-', label='Biofilm Reduction % (vs Initial)')
        axs[1].axhline(target_reduction_min, color='grey', linestyle='--', label=f'{target_reduction_min}% Target')
        axs[1].axhline(target_reduction_max, color='grey', linestyle='--')
        axs[1].set_ylabel("Reduction (%)")
        axs[1].set_title("Biofilm Reduction Over Time")
        axs[1].set_ylim(bottom=max(-10, stats_df['biofilm_pct_reduction_vs_initial'].min() - 10), top=110)
        axs[1].legend()
        axs[1].grid(True, linestyle=':')

        # Plot 3: Probiotic States (Stacked Area)
        probiotic_state_cols = ['probiotic_live_lime', 'probiotic_live_yellow',
                                'probiotic_live_cyan', 'probiotic_live_magenta', 'probiotic_dead']
        # MODIFIED: Changed labels to reflect functional states
        probiotic_state_labels = ['Killing (DNase+IL10)', 'Killing (DNase only)', 'Secreting IL10 only', 'Idle', 'Dead']
        probiotic_state_colors = ['lime', 'yellow', 'cyan', 'magenta', 'black'] # Colors remain the same
        # Ensure columns exist before plotting
        probiotic_data_to_plot = stats_df[probiotic_state_cols].copy()

        axs[2].stackplot(stats_df['time_min'], probiotic_data_to_plot.T, labels=probiotic_state_labels, colors=probiotic_state_colors, alpha=0.8)
        axs[2].set_ylabel("Probiotic Cell Count")
        axs[2].set_title("Probiotic Cell States Over Time")
        axs[2].legend(loc='upper left')
        axs[2].set_xlabel("Time (min)")
        axs[2].grid(True, linestyle=':')
        axs[2].set_ylim(bottom=0)

        # Adjust layout and save
        fig.tight_layout()
        plot_full_path = os.path.join(output_directory, output_plot_file) # Save inside output dir
        plt.savefig(plot_full_path, dpi=150)
        plt.close(fig)
        print(f"Successfully saved summary plots to: {plot_full_path}")

    except Exception as e:
        print(f"\nError generating plots: {e}")
        print("Ensure matplotlib is installed.")

else:
    print("\nNo data processed, CSV file and plots not saved.")

# --- Print Final Summary ---
if all_stats_data:
    final_biofilm_total = final_biofilm_live + final_biofilm_dead
    if final_biofilm_total > 0:
        percent_biofilm_live_final = (final_biofilm_live / final_biofilm_total) * 100.0
        percent_biofilm_dead_final = (final_biofilm_dead / final_biofilm_total) * 100.0
    else:
        percent_biofilm_live_final = 0.0 if final_biofilm_live == 0 else np.nan
        percent_biofilm_dead_final = 0.0 if final_biofilm_dead == 0 else np.nan

    print("\n--- Simulation Final Statistics ---")
    print(f"Initial Biofilm Cells:      {int(initial_biofilm_count)}")
    print(f"Initial Probiotic Cells:    {int(initial_probiotic_count)}")
    print(f"Initial Total Cells:        {int(initial_biofilm_count + initial_probiotic_count)}")
    print("-" * 30)
    print(f"Final Biofilm Cells (Live): {final_biofilm_live}")
    print(f"Final Biofilm Cells (Dead): {final_biofilm_dead}")
    print(f"  (Total Final Biofilm:    {final_biofilm_total})")
    print("-" * 30)
    # Updated print labels for probiotic states
    print(f"Final Probiotic (Killing, DNase+IL10): {final_probiotic_lime}")
    print(f"Final Probiotic (Killing, DNase only): {final_probiotic_yellow}")
    print(f"Final Probiotic (IL10 only):         {final_probiotic_cyan}")
    print(f"Final Probiotic (Idle):              {final_probiotic_magenta}")
    print(f"Final Probiotic (Dead):              {final_probiotic_dead}")
    probiotic_live_sum = final_probiotic_lime + final_probiotic_yellow + final_probiotic_cyan + final_probiotic_magenta
    print(f"  (Sum Probiotic Live: {probiotic_live_sum})")
    print("-" * 30)
    print(f"Final Total Cells:          {final_total_cells}")
    print("-" * 30)
    if not np.isnan(percent_biofilm_dead_final):
        print(f"Biofilm Live Pct (Live/Total Final): {percent_biofilm_live_final:.2f}%")
        print(f"Biofilm Reduction (Dead/Total Final):{percent_biofilm_dead_final:.2f}%")
    else:
         print("Biofilm percentages undefined (no biofilm cells found in final state).")
    if initial_biofilm_count > 0:
         percent_biofilm_killed_vs_initial = 100.0 * (initial_biofilm_count - final_biofilm_live) / initial_biofilm_count
         print(f"Biofilm Reduction (vs Initial):    {percent_biofilm_killed_vs_initial:.2f}%")

    print("-" * 30)
print("Done.")
