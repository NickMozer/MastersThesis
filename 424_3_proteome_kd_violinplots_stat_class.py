# Violin plots of all 281 knockdown proteomes  with Mann Whitney P in comparison to the whole proteome dataset

# Import libraries
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import numpy as np

# Initialize date and output path
input_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(input_dir)
output_dir_classifications = input_dir

# Load fc, rsd tables
fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
kd_columns=[col for col in fc_table.columns if "geneName" not in col]
fc_values_general=fc_table.iloc[:, 1:].values.flatten()
distribution_classifications = {}

def format_p_value(p_value):
    """Show first two or three non-zero decimals."""
    if np.isnan(p_value):
        return "nan"
    if p_value > 0.05:
        p_wert_dezimal="%.2f"% p_value
    elif p_value < 0.01:
        p_wert_dezimal=rf'${format(p_value, ".2e").split("e")[0]} \times 10^{{{int(format(p_value, ".2e").split("e")[1])}}}$'
    elif p_value <= 0.05:
        p_wert_dezimal="%.3f"% p_value  
    return p_wert_dezimal  
for kd in kd_columns:

    # Reshape the data for boxplot
    distribution_classifications[kd] = ("Insufficient Data", "NoType", "NoMedian", "Noq1", "Noq3", "Nomin", "Nomax", "Nohighfcs", "Nolowfcs")
    fc_values_kd = fc_table.loc[:, kd].values.flatten()

    # Mann-Whitney-Test
    stat, p = mannwhitneyu(fc_values_general, fc_values_kd, alternative='two-sided')
    p_value=format_p_value(p)

    # Find number of fcs above 1.0 or below -1.0
    greater_than_1_0 = (fc_values_kd > 1.0).sum(axis=0)
    lower_than_1_0 = (fc_values_kd < -1.0).sum(axis=0)

    fig = plt.figure(figsize=(5, 10)) # Adjusted figure size for potentially more metabolites
    ax = fig.add_subplot(1, 1, 1)

    # Violin width relative to figure width
    violin_width = 0.8
    # Make explicit x-ticks and limits for enough space
    ax.set_xticks([0])
    ax.set_xlim(-0.75, 0.75) # Space left and right

    # Violin Plot for non-pathway related responses (gray)
    sns.violinplot(data=fc_values_kd,
                    color='darkgray', alpha=0.5, linewidth=0.5, zorder=1, width=violin_width, ax=ax)
    
    # Get boxplot statistics
    median = fc_table[kd].median()
    q1 = fc_table[kd].quantile(0.25)
    q3 = fc_table[kd].quantile(0.75)
    whisker_low = fc_table[kd].min()
    whisker_high = fc_table[kd].max()

    # Add annotations
    plt.text(0.3, median, f'Median: {median:.2f}', verticalalignment='center', color='blue')
    plt.text(0.3, q1, f'25%: {q1:.2f}', verticalalignment='center', color='dodgerblue')
    plt.text(0.3, q3, f'75%: {q3:.2f}', verticalalignment='center', color='dodgerblue')
    plt.text(0.3, whisker_low, f'Min: {whisker_low:.2f}', verticalalignment='center', color='cornflowerblue')
    plt.text(0.3, whisker_high, f'Max: {whisker_high:.2f}', verticalalignment='center', color='cornflowerblue')

    plt.xlabel(f'{kd}\n(n={len(fc_values_kd)}; p={p_value})')
    plt.ylabel('Log2 fold change')

    plt.title("")
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    plt.tick_params(axis='y', which='both', right=False, left=True, labelleft=True)
    for pos in ['right', 'top']:
        plt.gca().spines[pos].set_visible(False)
    plt.xticks(ha='right')

    plt.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'

    # Save in two folders to prefilter more relevant pathways
    if p > 0.05:
        output_dir=os.path.join(input_dir, f"./kd_violins_p/not/27_proteome_{kd}_violin.png")
        distribution_classifications[kd] = (p_value, "ns", median, q1, q3, whisker_low, whisker_high, greater_than_1_0, lower_than_1_0)
    elif p < 0.00001:
        output_dir=os.path.join(input_dir, f"./kd_violins_p/sig5/27_proteome_{kd}_violin.png")
        distribution_classifications[kd] = (p_value, "*****", median, q1, q3, whisker_low, whisker_high, greater_than_1_0, lower_than_1_0)
    elif p < 0.0001:
        output_dir=os.path.join(input_dir, f"./kd_violins_p/sig4/27_proteome_{kd}_violin.png")
        distribution_classifications[kd] = (p_value, "****", median, q1, q3, whisker_low, whisker_high, greater_than_1_0, lower_than_1_0)
    elif p < 0.001:
        output_dir=os.path.join(input_dir, f"./kd_violins_p/sig3/27_proteome_{kd}_violin.png")
        distribution_classifications[kd] = (p_value, "***", median, q1, q3, whisker_low, whisker_high, greater_than_1_0, lower_than_1_0)
    elif p < 0.01:
        output_dir=os.path.join(input_dir, f"./kd_violins_p/sig2/27_proteome_{kd}_violin.png")
        distribution_classifications[kd] = (p_value, "**", median, q1, q3, whisker_low, whisker_high, greater_than_1_0, lower_than_1_0)
    elif p < 0.05:
        output_dir=os.path.join(input_dir, f"./kd_violins_p/sig1/27_proteome_{kd}_violin.png")
        distribution_classifications[kd] = (p_value, "*", median, q1, q3, whisker_low, whisker_high, greater_than_1_0, lower_than_1_0)
    plt.savefig(output_dir)
    plt.close()

    print(f"Violin plots for {kd} saved to {output_dir}")

# Save all distribution classifications to a xlsx file
classification_df = pd.DataFrame.from_dict(distribution_classifications, orient='index', columns=['P_value', 'Significance', 'Median', 'Lower_quartile', 'Higher_quartile', 'Minimum', 'Maximum', 'Num_fcs>1', 'Num_fcs<-1'])
classification_df.index.name = 'Knockdown'
output_path_xlsx = os.path.join(output_dir_classifications, "knockdown_significances.xlsx")
classification_df.to_excel(output_path_xlsx)