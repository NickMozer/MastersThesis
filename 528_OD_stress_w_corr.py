# Make a scatterplot of the OD(t=6.5h) per numbers of foldchange in- and decreases of stress proteins per knockdown

# Import Libraries
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
import multiprocessing

num_cores = multiprocessing.cpu_count()
jobs_n=num_cores-1
print(f"Number of free CPU cores: {num_cores}, use {jobs_n}")

# Initialize output path
output_excel_path = "09_stress_responses.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)

# Load fc, regulon tables
fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
# put genes from RDBECOLISFC00007_SigmulonGenes in Tabelle1
stress_table=pd.read_excel("./Stress_proteins_list.xlsx", "Tabelle1")

knockdown_genes = [col for col in fc_table.columns if col != 'geneName']

fc_table_2=fc_table.set_index("geneName")
stress_set=set(stress_table["stress_proteins"])
stress_in_index_list=list(stress_set.intersection(fc_table["geneName"]))
stress_df=fc_table_2.loc[stress_in_index_list, knockdown_genes]

# Save to Excel
with pd.ExcelWriter(output_excel_path) as writer:
    stress_df.to_excel(writer, sheet_name="stress_responses")

cutoff = 1.0
kd_nr = 297
re_nr = 1978
s1 = 64 # dot size
f_size = 30 # font size

# Initialize output path
output_excel_path = f"20_high_fc_kds_OD_{cutoff}.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_path = os.path.join(output_dir, output_dir, output_excel_path)

# Load dataframes
df_ODs = pd.read_excel("./TableS1-strains-info-growth.xlsx", sheet_name="Table S1", index_col=0)
df_ODs = df_ODs.rename(columns={'Target': 'Knockdown'})

df_stress = pd.read_excel("./09_stress_responses.xlsx", sheet_name="stress_responses")

# Identify the columns representing the knockdowns
kd_columns = df_stress.columns[df_stress.columns != 'gene']
stress_data_only = df_stress[kd_columns]

# Ensure the data is numeric.
stress_data_only = stress_data_only.apply(pd.to_numeric, errors='coerce')

# Count the number of stress responses for each knockdown based on the cutoff.
# This is done using boolean masking (e.g., stress_data_only > cutoff) followed by .sum() along the rows (axis=0), which counts the True values.

# 1. Count positive fold changes (> cutoff)
pos_counts = (stress_data_only > cutoff).sum(axis=0)
pos_counts_df = pos_counts.reset_index()
pos_counts_df.columns = ['Knockdown', f'above{cutoff}']

# 2. Count negative fold changes (< -cutoff)
neg_counts = (stress_data_only < -cutoff).sum(axis=0)
neg_counts_df = neg_counts.reset_index()
neg_counts_df.columns = ['Knockdown', f'below-{cutoff}']

# 3. Count absolute fold changes (abs() > cutoff)
abs_counts = (stress_data_only.abs() > cutoff).sum(axis=0)
abs_counts_df = abs_counts.reset_index()
abs_counts_df.columns = ['Knockdown', f'both{cutoff}']

# Merge the new count dataframes with the ODs dataframe.
kd_table_pos = pd.merge(pos_counts_df, df_ODs[['Knockdown', 'OD(t=6.5h)']], on='Knockdown')
kd_table_neg = pd.merge(neg_counts_df, df_ODs[['Knockdown', 'OD(t=6.5h)']], on='Knockdown')
kd_table_abs = pd.merge(abs_counts_df, df_ODs[['Knockdown', 'OD(t=6.5h)']], on='Knockdown')

# Scatter plot of kd_table
scatter_tables = [kd_table_pos, kd_table_neg, kd_table_abs]
for i, table in enumerate(scatter_tables):
    # Skip if the table is empty
    if table.empty:
        print(f"Warning: Table {i} is empty, skipping plotting.")
        continue

    if i == 0:
        table_name = "OD_no_stress_fc_pos"
        x_col = f'above{cutoff}'
        x_label = f"Stress proteins with log2 FC > {cutoff}"
    elif i == 1:
        table_name = "OD_no_stress_fc_neg"
        x_col = f'below-{cutoff}'
        x_label = f"Stress proteins with log2 FC < -{cutoff}"
    else:
        table_name = "OD_no_stress_fc_abs"
        x_col = f'both{cutoff}'
        x_label = f"NStress proteins with log2 FC > {cutoff}"

    x = table[x_col]
    y = table["OD(t=6.5h)"]
    
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(x, y,
               s=s1,
               marker='o',
               edgecolor='black',
               facecolor='black',
               alpha=1.0)
    
    # Calculate maximum values for both axes
    max_x = np.nanmax(x)
    max_y = np.nanmax(y)
    
    global_min = 0

    # Add correlation line
    # Perform linear regression
    mask = ~np.isnan(x) & ~np.isnan(y)
    if mask.any(): # Only perform regression if there are valid data points
        slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask], y[mask])
        
        # Create the regression line
        line_x = np.array([global_min, max_x])
        line_y = slope * line_x + intercept
        
        ax.plot(line_x, line_y, color='k', linestyle='-', linewidth=1, label=f'R² = {r_value**2:.2f}')
        ax.legend(fontsize=f_size) # Display the R-squared value as a legend

    # Determine ax ranges
    ax.set_xlim(global_min, max_x + 1) # Add a bit of padding
    ax.set_ylim(global_min, max_y + 0.1)

    ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labelsize=f_size)
    ax.tick_params(axis='y', which='both', right=False, left=True, labelleft=True, labelsize=f_size)
    for pos in ['right', 'top']:
        ax.spines[pos].set_visible(False)
    
    ax.xaxis.set_label_text(x_label, fontsize=f_size)
    ax.yaxis.set_label_text("OD(t=6.5h)", fontsize=f_size)
    
    plt.rcParams['svg.fonttype'] = 'none'
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"20_{table_name}_{cutoff}.svg"))
    plt.savefig(os.path.join(output_dir, f"20_{table_name}_{cutoff}.png"))
    plt.close(fig)

# Safe to Excel
with pd.ExcelWriter(output_path) as writer:
    kd_table_pos.to_excel(writer, sheet_name=f'Kd > {cutoff} FC_OD', index=False)
    kd_table_neg.to_excel(writer, sheet_name=f'Kd < -{cutoff} FC_OD', index=False)
    kd_table_abs.to_excel(writer, sheet_name=f'Abs_Kd > {cutoff} FC_OD', index=False)

print(f"Tables and figures saved to {output_dir}.")