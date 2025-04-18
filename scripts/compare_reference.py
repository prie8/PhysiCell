import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd # For organizing data

# Descriptive labels for your personalized smart bacteria simulation conditions
smart_bacteria_conditions = [
    "Smart Bacteria (Low HbA1c Sim)",
    "Smart Bacteria (Avg HbA1c Sim)",
    "Smart Bacteria (High HbA1c Sim)"
]
# Corresponding MEAN % biofilm MASS reduction at **120 hours**
smart_bacteria_reduction = [66.64, 62.12, 57.92]

# Vancomycin benchmark data (Ensure benchmark is also relevant for 120h or state timeframe)
vancomycin_label = "Vancomycin"
vancomycin_reduction = 30 # Example benchmark value
# --- End of data section ---


# Combine data for plotting
treatments = smart_bacteria_conditions + [vancomycin_label]
reduction_values = smart_bacteria_reduction + [vancomycin_reduction]

# Create a 'Type' column to distinguish treatments for coloring
treatment_types = ['Smart Bacteria'] * len(smart_bacteria_conditions) + ['Vancomycin']

# Create a Pandas DataFrame
df = pd.DataFrame({
    'Simulation Condition / Benchmark': treatments,
    'Biofilm Mass Reduction (%)': reduction_values,
    'Treatment Type': treatment_types
})

# --- Plotting ---
plt.figure(figsize=(10, 7))

# Create the bar plot using Seaborn
barplot = sns.barplot(
    x='Simulation Condition / Benchmark',
    y='Biofilm Mass Reduction (%)',
    data=df,
    hue='Treatment Type',
    palette=['#1f77b4', '#ff7f0e'], # Example: Blue for Smart Bacteria, Orange for Vanco
    dodge=False
    )

# Add labels and title (using updated metric name and time)
plt.xlabel("Simulation Condition / Benchmark", fontsize=12)
# *** UPDATED Y-AXIS LABEL BELOW ***
plt.ylabel("Predicted Biofilm Mass Reduction (%) at 120 hours", fontsize=12)
plt.title("Comparison of Predicted Biofilm Reduction Efficacy (120h)", fontsize=14) # Updated title too
plt.ylim(0, 100)

# Rotate x-axis labels if they overlap
plt.xticks(rotation=10, ha='right', fontsize=10)

# Add percentage values on top of bars
for p in barplot.patches:
    barplot.annotate(f"{p.get_height():.0f}%",
                   (p.get_x() + p.get_width() / 2., p.get_height()),
                   ha = 'center', va = 'center',
                   xytext = (0, 9),
                   textcoords = 'offset points',
                   fontsize=10, color='black')

# Adjust legend
plt.legend(title='Treatment Type', loc='upper left')

# Ensure layout is neat
plt.tight_layout()

# Save the figure
plt.savefig("comparison_biofilm_reduction_120h.png", dpi=300) # Updated filename

# Display the plot
plt.show()
