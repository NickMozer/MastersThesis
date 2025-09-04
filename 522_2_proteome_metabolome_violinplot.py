# Violin plots of proteome vs metabolome

# Import libraries
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Initialize date and output path
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)

# Load fc, rsd tables
fc_table_prot = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
fc_table_metab = pd.read_excel("./219_metabolome_fc_table_low_rsds.xlsx", "Sheet1")

# Reshape the data for violin plot
# Extract and flatten the numerical fold change values from both tables
fc_values_prot = fc_table_prot.iloc[:, 1:].values.flatten()
fc_values_metab = fc_table_metab.iloc[:, 1:].values.flatten()

print(f"Proteome FC values shape: {fc_values_prot.shape}")
print(f"Metabolome FC values shape: {fc_values_metab.shape}")

# Create a DataFrame in "long" format suitable for seaborn.violinplot
# Pad the shorter array with NaNs to make lengths equal for concatenation
max_len = max(len(fc_values_prot), len(fc_values_metab))

# Use pd.Series to ensure proper NaN padding and alignment
s_prot = pd.Series(fc_values_prot)
s_metab = pd.Series(fc_values_metab)

# Concatenate the series, they will be padded with NaN automatically if different lengths
results_df = pd.DataFrame({
    'response_fc': pd.concat([s_prot, s_metab], ignore_index=True),
    'position': ['prot'] * len(s_prot) + ['metab'] * len(s_metab)
})

print(results_df.head())
print(results_df.tail())

# Violinplot for prot vs. metab
plt.figure(figsize=(8, 10))
sns.violinplot(x='position', y='response_fc', data=results_df, order=['prot', 'metab'], color="lightgray")
plt.ylabel('Log2 fold change')
plt.xlabel('')

# Drop NA values before calculating length for labels
prot_fcs = results_df[results_df['position'] == 'prot']['response_fc'].dropna()
metab_fcs = results_df[results_df['position'] == 'metab']['response_fc'].dropna()
gewuenschte_ticks = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
plt.yticks(gewuenschte_ticks)
plt.xticks([0, 1], [f'Proteome\n(n={len(prot_fcs)})', f'Metabolome\n(n={len(metab_fcs)})'])
sns.despine(top=True, right=True)
plt.tight_layout()
plt.rcParams['svg.fonttype'] = 'none'
plt.savefig(os.path.join(output_dir, f"24_proteome_fcs_violin.png"))
plt.close()

# Get boxplot statistics
median = s_prot.median()
q1 = s_prot.quantile(0.25)
q3 = s_prot.quantile(0.75)
whisker_low = s_prot.min()
whisker_high = s_prot.max()
print(f"Proteome median: {median}, q1: {q1}, q3: {q3}, whisker_low: {whisker_low}, whisker-high: {whisker_high}")

# Get boxplot statistics
median = s_metab.median()
q1 = s_metab.quantile(0.25)
q3 = s_metab.quantile(0.75)
whisker_low = s_metab.min()
whisker_high = s_metab.max()
print(f"Metabolome median: {median}, q1: {q1}, q3: {q3}, whisker_low: {whisker_low}, whisker-high: {whisker_high}")
print(f"Violin plots for proteome and metabolome saved to {output_dir}")