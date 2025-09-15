# Violinplots of response genes up- or downstream in the same operon as a knockdown/target gene (mannwhitney test)
# Also plot all fcs in a violin plot for reference
# Find out medians of each group

# Import libraries
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu

# Initialize output path
output_excel_path = "032_operon_analysis.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_path = os.path.join(output_dir, output_excel_path)

# Load fc, operon tables
fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
operon_table=pd.read_excel("./operons_ecoli.xlsx", "Tabelle1")

# Filter for operons with more than one gene
filtered_operon_table = pd.DataFrame(columns=operon_table.columns)
for i, row in operon_table.iterrows():
    if any("-" in str(value) for value in row):
        filtered_operon_table = pd.concat([filtered_operon_table, row.to_frame().T], ignore_index=True)

# Knockdown/target columns
knockdown_columns = fc_table.columns.difference(['geneName'])

# Define functions
def format_p_value(p_value):
    if np.isnan(p_value):
        return "nan"
    if p_value > 0.05:
        p_wert_dezimal="%.2f"% p_value
    elif p_value < 0.01:
        p_wert_dezimal=rf'${format(p_value, ".2e").split("e")[0]} \times 10^{{{int(format(p_value, ".2e").split("e")[1])}}}$'
    elif p_value <= 0.05:
        p_wert_dezimal="%.3f"% p_value  
    return p_wert_dezimal  

def calculate_mannwhitney_p_value(group1, group2):
    if len(group1) < 2 or len(group2) < 2:
        return 1.0
    stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
    return p

def calculate_and_annotate_significance_single_plot(group1, group2):
    if len(group1) < 2 or len(group2) < 2:
        return "ns"

    stat, p = mannwhitneyu(group1, group2, alternative='two-sided')

    # Significance symbols
    if p >= 0.05:
        significance = "ns"
    elif p < 0.00001:
        significance = "*****"
    elif p < 0.0001:
        significance = "****"
    elif p < 0.001:
        significance = "***"
    elif p < 0.01:
        significance = "**"
    elif p < 0.05:
        significance = "*"
    return significance

# Analyse knockdowns in an operon, differ between upstream and downstream responses
results = []
for index, operon_row in filtered_operon_table.iterrows():
    operon_genes = operon_row['CsiR'].split('-')
    for knocked_down_gene in operon_genes:
        if knocked_down_gene in knockdown_columns:
            for gene in operon_genes:
                if gene != knocked_down_gene:
                    try:
                        knockdown_fc_values = fc_table.loc[fc_table['geneName'] == gene, knocked_down_gene].values[0]

                        position = "upstream" if operon_genes.index(gene) < operon_genes.index(knocked_down_gene) else "downstream"

                        results.append({
                            'knockdown_gene': knocked_down_gene,
                            'response_gene': gene,
                            'response_fc': knockdown_fc_values,
                            'position': position
                        })
                        print(f"Gene {gene} found in fc_table.")
                    except IndexError:
                        print(f"Gene {gene} not found.")  

# Append all the fc values for third violin plot, this takes the most time
response_genes=fc_table['geneName']
for knockdown in knockdown_columns:
    for response in response_genes:
        if response == knockdown:
            continue
        elif knockdown == "trpD" and response == "trpGD":
            continue
        knockdown_fc_values = fc_table.loc[fc_table['geneName'] == response, knockdown].values[0]
        results.append({
            'knockdown_gene': f"{knockdown}_01",
            'response_gene': f"{response}_01",
            'response_fc': knockdown_fc_values,
            'position': "Fold changes of non-targets"
        })
    print(f"Knockdown {knockdown} found for all fcs.")  

results_df = pd.DataFrame(results)

# Violinplot for Fold changes of non-targets, Upstream vs. Downstream with mannwhitney test
plt.figure(figsize=(10, 10))
sns.violinplot(x='position', y='response_fc', data=results_df, order=['Fold changes of non-targets', 'upstream', 'downstream'], color="lightgray")
plt.ylabel('Log2 fold change')
plt.xlabel('Position in operon relative to knockdown')

# mannwhitney-Test 
upstream_fcs = results_df[results_df['position'] == 'upstream']['response_fc'].dropna()
downstream_fcs = results_df[results_df['position'] == 'downstream']['response_fc'].dropna()
all_fcs= results_df[results_df['position'] == 'Fold changes of non-targets']['response_fc'].dropna()

upstream_median=upstream_fcs.median()
downstream_median=downstream_fcs.median()
all_median=all_fcs.median()
print(f"Upstream median: {upstream_median}, downstream median: {downstream_median}, all fcs median: {all_median}")

# Annotation of Signifikance bar
y_max = max(upstream_fcs.max() if not upstream_fcs.empty else -np.inf,
            downstream_fcs.max() if not downstream_fcs.empty else -np.inf,
            all_fcs.max() if not all_fcs.empty else -np.inf)
y_min = min(upstream_fcs.min() if not upstream_fcs.empty else np.inf,
            downstream_fcs.min() if not downstream_fcs.empty else np.inf,
            all_fcs.min() if not all_fcs.empty else np.inf)
y_range = y_max - y_min if y_max != -np.inf and y_min != np.inf else 1

if len(upstream_fcs) >= 2 and len(downstream_fcs) >= 2:
    mannwhitney_p = calculate_mannwhitney_p_value(upstream_fcs, downstream_fcs)
    significance = calculate_and_annotate_significance_single_plot(upstream_fcs, downstream_fcs)
    formatted_p = format_p_value(mannwhitney_p)

    y_bar = y_max + 0.05 * y_range
    y_text = y_bar + 0.02 * y_range
    x_positions = [1, 2]  # Positions of Violinplots (upstream=1, downstream=2)

    plt.plot(x_positions, [y_bar, y_bar], color='black', linewidth=1.5)
    plt.text(1.5, y_text, f"{significance} (P={formatted_p})", ha='center', va='bottom', color='black')

else:
    plt.text(0.5, y_max + 0.1 * y_range if y_max != -np.inf else 1, "Nicht genügend Daten für mannwhitney-Test", ha='center', color='red')

plt.xticks([0, 1, 2], [f'Fold changes of non-targets\n(n={len(all_fcs)})', f'Upstream\n(n={len(upstream_fcs)})', f'Downstream\n(n={len(downstream_fcs)})'])
sns.despine(top=True, right=True)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, f"032_Operon_violinplot_all_fcs_upstream_downstream_mannwhitney.png"))
plt.savefig(os.path.join(output_dir, f"032_Operon_violinplot_all_fcs_upstream_downstream_mannwhitney.svg"))
plt.close()

# Save table to Excel
with pd.ExcelWriter(output_path) as writer:
    filtered_operon_table.to_excel(writer, sheet_name=f"filtered_operons", index=False)
    results_df.to_excel(writer, sheet_name=f"operon_fcs", index=False)
print(f"Table and figure saved in {output_dir}.")
