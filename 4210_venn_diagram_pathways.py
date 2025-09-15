# Venn diagram of Ecocyc subpathways, response pathways and target pathways

# Import libraries
import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import os
import math

# Define output directory
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_excel_path = os.path.join(output_dir, f"05_pathways_venn.xlsx")
output_path = os.path.join(output_dir, output_excel_path)

# Load the pathway table
pathway_table = pd.read_excel("./20250327_pathways_reactants_ecoli.xlsx", "subpathways")
pathway_table = pathway_table.rename(columns={'GeneName': 'geneName', 'Pathway': 'pathway'})

# Load the fold change table
fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
knockdown_genes = [col for col in fc_table.columns if col != 'geneName']

# 1. Get the set of all unique subpathways from Ecocyc
ecocyc_subpathways = set(pathway_table['pathway'].unique())
ecocyc_subpathways = {
    pathway for pathway in ecocyc_subpathways if not isinstance(pathway, float) or not math.isnan(pathway)
}
print(f"Anzahl der eindeutigen Subpathways in Ecocyc: {len(ecocyc_subpathways)}")

# 2. Get the set of pathways associated with the response genes (all genes in fc_table)
def get_unique_genes(series):
    return list(set(series))

response_gene_pathways_df = pd.merge(fc_table[['geneName']], pathway_table[['geneName', 'pathway']], on='geneName', how='inner')
response_gene_pathways_with_genes = response_gene_pathways_df.groupby('pathway')['geneName'].apply(get_unique_genes).reset_index(name='genes')
response_gene_pathways = set(response_gene_pathways_df['pathway'].unique())
response_gene_pathways = {
    pathway for pathway in response_gene_pathways if not isinstance(pathway, float) or not math.isnan(pathway)
}
print(f"Anzahl der eindeutigen Pathways in den response-Genen: {len(response_gene_pathways)}")

# 3. Get pathways and genes for knockdown genes
knockdown_gene_pathway_gene_df = pathway_table[pathway_table['geneName'].isin(knockdown_genes)].groupby('pathway')['geneName'].apply(get_unique_genes).reset_index(name='genes')
knockdown_gene_pathways = set(knockdown_gene_pathway_gene_df['pathway'])
knockdown_gene_pathways = {
    pathway for pathway in knockdown_gene_pathways if not isinstance(pathway, float) or not math.isnan(pathway)
}
print(f"Anzahl der eindeutigen Pathways in den Knockdown-Genen: {len(knockdown_gene_pathways)}")

# Create the Venn diagram
plt.figure(figsize=(20, 10))
plt.rcParams['font.size'] = 24
colors = ("royalblue", "cornflowerblue", "cyan")
v = venn3([ecocyc_subpathways, response_gene_pathways, knockdown_gene_pathways],
          set_labels=('Ecocyc Subpathways', 'Pathways der response-Gene', 'Pathways der Knockdown-Gene'),
          set_colors=colors)

# Set alpha only if the patch exists
if v.get_patch_by_id('100') is not None:
    v.get_patch_by_id('100').set_alpha(0.7)
if v.get_patch_by_id('010') is not None:
    v.get_patch_by_id('010').set_alpha(0.7)
if v.get_patch_by_id('001') is not None:
    v.get_patch_by_id('001').set_alpha(0.7)

# Adjust color of label texts
if v.get_label_by_id('A'):
    v.get_label_by_id('A').set_text(f"Ecocyc subpathways ({len(ecocyc_subpathways)})")
    v.get_label_by_id('A').set_color(colors[0])
if v.get_label_by_id('B'):
    v.get_label_by_id('B').set_text(f"Response-gene pathways ({len(response_gene_pathways)})")
    v.get_label_by_id('B').set_color(colors[1])
if v.get_label_by_id('C'):
    v.get_label_by_id('C').set_text(f"Target-gene pathways({len(knockdown_gene_pathways)})")
    v.get_label_by_id('C').set_color(colors[2])

plt.title("Venn-diagram of pathways in Ale proteome data")

# Save the Venn diagram
plt.rcParams['svg.fonttype'] = 'none'
plt.savefig(os.path.join(output_dir, "05_venn_diagram_pathways.png"))
plt.savefig(os.path.join(output_dir, "05_venn_diagram_pathways.svg"))
plt.close()

# Save to Excel
with pd.ExcelWriter(output_excel_path) as writer:
    pd.Series(list(ecocyc_subpathways)).to_excel(writer, sheet_name=f"Ecocyc_pathways", index=False)
    response_gene_pathways_with_genes.to_excel(writer, sheet_name=f"Response_pathways", index=False)
    knockdown_gene_pathway_gene_df.to_excel(writer, sheet_name=f"Knockdown_pathways", index=False)
print(f"Venn-Diagramm gespeichert unter: {output_path}")
