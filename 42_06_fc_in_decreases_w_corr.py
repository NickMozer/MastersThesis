# Make a scatterplot of the numbers of foldchange in- and decreases per knockdown
# Make a second plot for knockdowns per response genes
# Include correlation line and Pearson R

# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats

kd_nr=297
re_nr=1978
s1=64 # dot size
f_size=30 # font size

# Initalize output path
output_excel_path = "02_fc_in_and_decreases.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_path = os.path.join(output_dir, output_excel_path)

df_fc = pd.read_excel("./005_log2_fcs_final.xlsx", sheet_name="Successfull_knockdowns", index_col=0)

# Finding the TOP Knockdowns with the most log2 fold change values > 1.0
greater_than_1_0 = (df_fc > 1.0).sum(axis=0).nlargest(kd_nr)
top_knockdowns_high_fc = pd.DataFrame({'Knockdown': greater_than_1_0.index, 'Log2 fold change > 1.0 count': greater_than_1_0.values})

# Finding the TOP Knockdowns with the most log2 fold change values < -1.0
lower_than_1_0 = (df_fc < -1.0).sum(axis=0).nlargest(kd_nr)
top_knockdowns_low_fc = pd.DataFrame({'Knockdown': lower_than_1_0.index, 'Log2 fold change < -1.0 count': lower_than_1_0.values})

kd_table = pd.merge(top_knockdowns_high_fc, top_knockdowns_low_fc, on='Knockdown')

# Finding the TOP Proteins (responses) with the most log2 fold change values > 1.0
greater_than_1_0_rows = (df_fc > 1.0).sum(axis=1).nlargest(re_nr)
top_responses_high_fc = pd.DataFrame({'Response': greater_than_1_0_rows.index, 'Log2 fold change > 1.0 count': greater_than_1_0_rows.values})

# Finding the TOP Proteins (responses) with the most log2 fold change values < -1.0
lower_than_1_0_rows = (df_fc < -1.0).sum(axis=1).nlargest(re_nr)
top_responses_low_fc = pd.DataFrame({'Response': lower_than_1_0_rows.index, 'Log2 fold change < -1.0 count': lower_than_1_0_rows.values})

resp_table = pd.merge(top_responses_high_fc, top_responses_low_fc, on='Response')

# Scatter plot of kd_table
scatter_tables = [kd_table, resp_table]
for i, table in enumerate(scatter_tables):
    x = table[f"Log2 fold change < -1.0 count"]
    y = table[f"Log2 fold change > 1.0 count"]
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

    global_max = max(max_x, max_y)*1.01

    global_min = 0

    # Add correlation line
    mask = ~np.isnan(x) & ~np.isnan(y)
    if mask.any(): # Only perform regression if there are valid data points
        slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask], y[mask])
        
        # Create the regression line
        # Use np.linspace to get points across the desired range for the line
        line_x = np.array([global_min, global_max])
        line_y = slope * line_x + intercept
        
        ax.plot(line_x, line_y, color='k', linestyle='-', linewidth=1, label=f'R² = {r_value**2:.2f}')
        ax.text(global_max * 0.5/16, global_max *14.5/16, f'R² = {r_value**2:.2f}', fontsize=f_size)

    if i == 0:
        table_name = "Number of high fold changes in response proteins per knockdown"
        table_name1 = "A"
    else:
        table_name = "Number of knockdowns with high fold change in response protein"
        table_name1 = "B"

    # Get the linewidth of one of the axes spines
    axis_linewidth = ax.spines['bottom'].get_linewidth()

    # Determine ax ranges
    ax.set_xlim(global_min, global_max)
    ax.set_ylim(global_min, global_max)

    ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
    ax.tick_params(axis='y', which='both', right=False, left=True, labelleft=True)
    for pos in ['right', 'top']:
        ax.spines[pos].set_visible(False)

    # Set unique IDs for the x and y axes
    ax.xaxis.set_label_text("Number of log2 fold changes < -1.0", fontsize=f_size)
    ax.yaxis.set_label_text("Number of log2 fold changes > 1.0", fontsize=f_size)
    ax.xaxis.label.set_gid("x_axis_" + table_name)
    ax.yaxis.label.set_gid("y_axis_" + table_name)
    ax.xaxis.set_ticklabels(ax.get_xticklabels(), gid="x_axis_ticks_" + table_name, fontsize=f_size)
    ax.yaxis.set_ticklabels(ax.get_yticklabels(), gid="y_axis_ticks_" + table_name, fontsize=f_size)
    ax.spines['bottom'].set_gid("x_axis_line_" + table_name)
    ax.spines['left'].set_gid("y_axis_line_" + table_name)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"02_{table_name}_corr.svg"))
    plt.savefig(os.path.join(output_dir, f"02_{table_name}_corr.png"))
    plt.close(fig)

# Safe to Excel
with pd.ExcelWriter(output_path) as writer:
    top_knockdowns_high_fc.to_excel(writer, sheet_name=f' Knockdowns > 1.0 FC')
    top_knockdowns_low_fc.to_excel(writer, sheet_name=f'Knockdowns < -1.0 FC')
    top_responses_high_fc.to_excel(writer, sheet_name=f'Responses > 1.0 FC')
    top_responses_low_fc.to_excel(writer, sheet_name=f'Responses < -1.0 FC')

print(f"Tables and figures saved to {output_dir}.")
