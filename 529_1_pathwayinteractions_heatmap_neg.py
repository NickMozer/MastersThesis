# Create heatmap of pathway responses with fraction of responding genes with log2-fcs below a cutoff of log2 fc -1.0 upon knockdown of genes
# Same gene as knockdown gene excluded for pathway response
# Colorbar height adjusted for heatmap with numbers for x labels (not possible for normal labels)
# Filter out pathways that have less than 1% response fractions in their highest response or knockdown pathway
# Filter out pathways that have less than 5% response fractions in their mean response or knockdown pathways
# Sort pathways by subsystem
# Clipping at lower values than 1.0 response possible

# Import libraries
import pandas as pd
import os
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec

num_cores_avail = multiprocessing.cpu_count()
num_cores=num_cores_avail-1
print(f"Number of free CPU cores: {num_cores_avail}, use {num_cores}")

# Log2 fold change cutoff
cutoff = -1.0
# Minimum mean per pathway row or column in heatmap
co2 = 0.05
# Minimum Max-value per pathway row or column in heatmap
mini = 0.01
clip_up=0.5

# Initialize output path
output_excel_path = "04_heatmap_data.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_path = os.path.join(output_dir, output_excel_path)

fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
pathway_table = pd.read_excel("./20250327_pathways_reactants_ecoli.xlsx", "subpathways")
pathway_table = pathway_table.rename(columns={'GeneName': 'geneName', 'Pathway': 'pathway'})
subsys_order_df=pd.read_excel("./Subsys_order.xlsx", "subsys_sorting")
merged_table = pd.merge(fc_table, pathway_table[['geneName', 'pathway']], on='geneName', how='inner')
knockdown_genes = [col for col in fc_table.columns if col != 'geneName']

# Load Operon-Table
operon_table = pd.read_excel("./operons_ecoli.xlsx", "Tabelle1")

# Filter the operon table to only include rows where 'CsiR' contains '-'
filtered_operon_table = operon_table[operon_table['CsiR'].astype(str).str.contains('-')].copy()

# function to test wether two genes lie within one operon
def are_genes_in_same_operon(gene1, gene2, operon_df):
    for _, row in operon_df.iterrows():
        if isinstance(row['CsiR'], str) and '-' in row['CsiR']:
            operon_genes = row['CsiR'].split('-')
            if gene1 in operon_genes and gene2 in operon_genes:
                return True
    return False

def get_unique_genes(series):
    return list(set(series))

knockdown_gene_pathway_gene_df = pathway_table[pathway_table['geneName'].isin(knockdown_genes)].groupby('pathway')['geneName'].apply(get_unique_genes).reset_index(name='kd_genes')

knockdown_gene_pathway_gene_df=knockdown_gene_pathway_gene_df.sort_values("kd_genes").reset_index(drop=True)

response_gene_pathways_df = pd.merge(fc_table[['geneName']], pathway_table[['geneName', 'pathway']], on='geneName', how='inner')
response_gene_pathways_with_genes = response_gene_pathways_df.groupby('pathway')['geneName'].apply(get_unique_genes).reset_index(name='resp_genes')
response_gene_pathways_with_genes = response_gene_pathways_with_genes.sort_values("resp_genes").reset_index(drop=True)

merged_kr_df= pd.merge(knockdown_gene_pathway_gene_df, response_gene_pathways_with_genes, on='pathway', how='outer')
merged_kr_df=merged_kr_df.sort_values("kd_genes").reset_index(drop=True)

def filter_and_merge_pathways(mergdf, kddf):
    rows_to_drop = []
    rows_to_add = []
    rows_to_add1 = []
    rows_to_add2 = []
    processed_indices = set()

    for index1, row1 in mergdf.iterrows():
        if index1 in processed_indices:
            continue

        kd_gene_list1 = row1['kd_genes']
        resp_gene_list1 = row1['resp_genes']
        pathway1 = row1['pathway']

        for index2, row2 in mergdf.iterrows():
            if index2 == index1:
                continue
            kd_gene_list2 = row2['kd_genes']
            pathway2 = row2['pathway']
            resp_gene_list2 = row2['resp_genes']

            if kd_gene_list1 == kd_gene_list2 and resp_gene_list1 == resp_gene_list2:
                # Summarize rows
                merged_pathway = f"{pathway1};{pathway2}"
                rows_to_add.append({'pathway': merged_pathway, 'kd_genes': kd_gene_list1, 'resp_genes': resp_gene_list1})
                rows_to_add1.append({'pathway': merged_pathway, 'kd_genes': kd_gene_list1})
                rows_to_add2.append({'pathway': pathway1, 'kd_genes': kd_gene_list1})
                rows_to_drop.extend([index1, index2])
                processed_indices.add(index2)
                break

    # Delete old rows, add new ones 
    mergdf_filtered = mergdf.drop(rows_to_drop).reset_index(drop=True)
    mergdf_final = pd.concat([mergdf_filtered, pd.DataFrame(rows_to_add)], ignore_index=True)
    kddf_filtered = kddf.drop(rows_to_drop).reset_index(drop=True)
    kddf_final_m = pd.concat([kddf_filtered, pd.DataFrame(rows_to_add1)], ignore_index=True)
    kddf_final_1 = pd.concat([kddf_filtered, pd.DataFrame(rows_to_add2)], ignore_index=True)
    return mergdf_final, kddf_final_m, kddf_final_1

merged_kr_filtered, knockdown_gene_pathway_gene_df_m_filtered, knockdown_gene_pathway_gene_df_1_filtered = filter_and_merge_pathways(merged_kr_df.copy(), knockdown_gene_pathway_gene_df.copy())
kd_pathways = set(knockdown_gene_pathway_gene_df_1_filtered['pathway'])
print(knockdown_gene_pathway_gene_df_1_filtered)

# Mapping knockdown-gene -> knockdown-pathway
kd_pathway_to_genes = knockdown_gene_pathway_gene_df_1_filtered.set_index('pathway')['kd_genes'].to_dict()
kd_mapping = pathway_table.set_index('geneName')['pathway'].to_dict()

kd_genes_in_pathways=list()
kd_genes_not_in_pathways=list()
for gene in knockdown_genes:
    if gene in kd_mapping:
        kd_genes_in_pathways.append(gene)
    else:
        kd_genes_not_in_pathways.append(gene)

print(f"Loaded all tables successfully, made dictionary of knockdown genes per pathway, {len(kd_genes_in_pathways)}, {len(kd_genes_not_in_pathways)}")
pathway_legend_results = {}
# Function to analyze a single pathway
def analyze_pathway(pathway_name, pathway_df, cutoff, fc_table):
    pathway_gene_names = pathway_df['geneName'].tolist()
    response = {}

    for kd_pathway in kd_pathways:
        num_all_relevant_genes=0
        num_all_responding_genes=0
        kd_pathway_gene_names = kd_pathway_to_genes[kd_pathway]
        for kd_gene in kd_pathway_gene_names:
            if kd_gene in fc_table.columns:
                # Filter genes that are in the fc table and not the knockdown gene
                relevant_pathway_genes = [
                    gene for gene in pathway_gene_names
                    if gene in fc_table['geneName'].values and gene != kd_gene 
                ]
                num_relevant_genes = len(relevant_pathway_genes)
                num_all_relevant_genes=num_all_relevant_genes+num_relevant_genes
                
                if num_relevant_genes > 0:
                    if kd_gene in filtered_operon_table:
                        # Count genes in relevant pathway with foldchanges below cutoff
                        relevant1_pathway_genes = [
                            gene for gene in pathway_gene_names
                            if gene in fc_table['geneName'].values and gene != kd_gene and not are_genes_in_same_operon(kd_gene, gene, filtered_operon_table)
                        ]
                        num_relevant1_genes = len(relevant_pathway_genes)
                        responding_genes = [
                            fc_table.loc[fc_table['geneName'] == gene, kd_gene].values[0]
                            for gene in relevant1_pathway_genes
                            if fc_table.loc[fc_table['geneName'] == gene, kd_gene].values[0] < cutoff
                        ]
                        num_responding_genes = len(responding_genes)
                        fraction_responding = num_responding_genes / num_relevant1_genes
                        response[kd_gene] = fraction_responding
                    else:
                        # Count genes in relevant pathway with foldchanges above cutoff
                        responding_genes = [
                            fc_table.loc[fc_table['geneName'] == gene, kd_gene].values[0]
                            for gene in relevant_pathway_genes
                            if fc_table.loc[fc_table['geneName'] == gene, kd_gene].values[0] < cutoff
                        ]
                        num_responding_genes = len(responding_genes)
                        num_all_responding_genes=num_all_responding_genes+num_responding_genes
                    
        fraction_responding = num_all_responding_genes / num_all_relevant_genes if num_all_relevant_genes > 0 else np.nan
        response[kd_pathway] = fraction_responding

    print(f"Pathway {pathway_name} for {cutoff} is done.")
    return pathway_name, response

# Analyze Pathways
results = Parallel(n_jobs=num_cores)(
    delayed(analyze_pathway)(name, df, cutoff, fc_table)
    for name, df in merged_table.groupby('pathway')
)
print("Interactions analysis done")
pathway_responses = dict(results)
print("Made heatmap data dict")

# Make heatmap data df
heatmap_data1 = pd.DataFrame.from_dict(pathway_responses, orient='index').reindex(columns=kd_pathways)
heatmap_data1.insert(0, 'pathway', heatmap_data1.index)
heatmap_data_excel_1 = heatmap_data1.copy()
heatmap_data_excel_1.to_excel("04_pathway_interactions_cutoff_-1.0.xlsx", sheet_name=f"fraction_resp_below_{cutoff}", index=False)
heatmap_data = pd.read_excel("04_pathway_interactions_cutoff_-1.0.xlsx", f"fraction_resp_below_{cutoff}", index_col="pathway")
print("Made and saved heatmap data df")

# Keep only pathways that were also knocked down
common_pathways_all = [p for p in kd_pathways if p in heatmap_data.index]
common_pathways=list(set(common_pathways_all))

# Heatmap data subset
heatmap_data_matched = heatmap_data.loc[common_pathways, common_pathways]
print("Matched heatmap data")

# Mean and max responses per row and column
row_means = heatmap_data_matched.mean(axis=1)
row_max = heatmap_data_matched.max(axis=1)
col_means = heatmap_data_matched.mean(axis=0)
col_max = heatmap_data_matched.max(axis=0)

# Identify pathway rows to keep
rows_to_keep_mean = row_means[row_means >= co2].index.tolist()
rows_to_keep_max = row_max[row_max > mini].index.tolist()
rows_to_keep = list(set(rows_to_keep_mean) | set(rows_to_keep_max))  # Unionize both conditions

# Identify pathway columns to keep
cols_to_keep_mean = col_means[col_means >= co2].index.tolist()
cols_to_keep_max = col_max[col_max > mini].index.tolist()
cols_to_keep = list(set(cols_to_keep_mean) | set(cols_to_keep_max))  # Unionize both conditions

# Find common pathways to keep in rows and columns
common_kept_pathways_01 = list(set(rows_to_keep) | set(cols_to_keep))
print(len(common_kept_pathways_01))

# Heatmap data subset
heatmap_data_matched_filtered_01 = heatmap_data.loc[common_kept_pathways_01, common_kept_pathways_01]
print("Matched heatmap data filtered")

# Sort Pathway order by Subsystem
merged_pathways=pd.merge(heatmap_data_matched_filtered_01, pathway_table[['pathway', 'subSys']], on='pathway', how="inner")
df_without_gene_not_found = merged_pathways[~merged_pathways.apply(lambda row: row.astype(str).str.contains('Gene Not found').any(), axis=1)]
unique_pathways_df = df_without_gene_not_found.groupby('pathway')['subSys'].first().reset_index()
unique_pathways_df_sorted=pd.merge(subsys_order_df[['sys_number', 'subSys']], unique_pathways_df[['pathway', 'subSys']], on="subSys", how="inner")
unique_pathways_df_sorted =unique_pathways_df_sorted.sort_values(by="sys_number").reset_index(drop=True)
common_kept_pathways=list(unique_pathways_df_sorted["pathway"])
heatmap_data_matched_filtered = heatmap_data.loc[common_kept_pathways, common_kept_pathways]
print("Matched heatmap data ordered by Subsystem")

# Fraction values above 0.3 set to 0.3
heatmap_data_clipped = heatmap_data_matched_filtered.clip(upper=clip_up)
print("Clipped heatmap data")

# Prepare for Excel
heatmap_data_excel_1 = heatmap_data.copy()
heatmap_data_excel_1.insert(0, 'pathway', heatmap_data_excel_1.index)
heatmap_data_excel_2 = heatmap_data_matched.copy()
heatmap_data_excel_2.insert(0, 'pathway', heatmap_data_excel_2.index)
heatmap_data_excel_3 = heatmap_data_matched_filtered.copy()
heatmap_data_excel_3.insert(0, 'pathway', heatmap_data_excel_3.index)
heatmap_data_excel_4 = heatmap_data_clipped.copy()
heatmap_data_excel_4.insert(0, 'pathway', heatmap_data_excel_4.index)

kd_genes_in_p_df=pd.DataFrame()
kd_genes_in_p_df["Knockdowns in Pathways"]=pd.DataFrame(kd_genes_in_pathways)
kd_genes_in_p_df["Knockdown Pathways"]=pd.DataFrame(kd_pathways)
kd_genes_in_p_df["Knockdowns not in Pathways"]=pd.DataFrame(kd_genes_not_in_pathways)

# Save to Excel
with pd.ExcelWriter(output_excel_path) as writer: 
    kd_genes_in_p_df.to_excel(writer, sheet_name="Knockdowns in Pathways", index=False)
    heatmap_data_excel_1.to_excel(writer, sheet_name=f"fraction_responses_over_{cutoff}", index=False)  
    heatmap_data_excel_2.to_excel(writer, sheet_name=f"pathway_interactions_matched", index=False)
    heatmap_data_excel_3.to_excel(writer, sheet_name=f"pathway_interactions_filtered", index=False)
    heatmap_data_excel_4.to_excel(writer, sheet_name=f"pathway_interactions_clip_{clip_up}", index=False)
print("Saved tables to Excel")

# Check if data to plot exists
if not heatmap_data_clipped.empty:
    # Heatmap Plot with gridspec for precise colorbar-placement
    plt.rcParams['font.size'] = 35/80*len(common_kept_pathways)
    fig = plt.figure(figsize=(max(10, heatmap_data_clipped.shape[1]), max(10, heatmap_data_clipped.shape[0])))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.05]) # Heatmap takes majority, colorbar smaller portion

    ax_heatmap = plt.subplot(gs[0])
    cbar_ax = plt.subplot(gs[1])

    sns.heatmap(heatmap_data_clipped, cmap='PiYG_r', center=0, annot=False, 
                linewidths=.1, linecolor="gray", cbar_kws={'label': 'Pathway response: fraction'}, square=True, ax=ax_heatmap, cbar_ax=cbar_ax)

    ax_heatmap.set_yticks(np.arange(len(common_kept_pathways)) + 0.5) # +0.5 centers labels
    ax_heatmap.set_yticklabels(common_kept_pathways, rotation=0) # Set labels explicitly
    ax_heatmap.set_xticks(np.arange(len(common_kept_pathways)) + 0.5)
    ax_heatmap.set_xticklabels(common_kept_pathways, rotation=90) # Set labels explicitly
    ax_heatmap.tick_params(axis='y', rotation=0)
    ax_heatmap.tick_params(axis='x', rotation=90)
    ax_heatmap.set_xlabel('Knockdown-Pathway')
    ax_heatmap.set_ylabel('Response-Pathway')

    plt.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(os.path.join(output_dir, f"04_pathway_interactions_heatmap_{cutoff}_and_{co2}_{mini}.png"))
    plt.savefig(os.path.join(output_dir, f"04_pathway_interactions_heatmap_{cutoff}_and_{co2}_{mini}.svg"))
    plt.close(fig)

else:
    print(f"No relevant data after filtering found. No heatmap made.")

unique_pathways_df_sorted.to_excel("04_ordered_heatmap_data_subsys.xlsx", "ordered_by_subsys")
pathway_number_mapping = {}
number_pathway_mapping = {}
number = 1

for pathway in common_kept_pathways:
    pathway_number_mapping[pathway] = number
    number_pathway_mapping[number] = pathway
    number += 1

df_legend = pd.DataFrame(list(number_pathway_mapping.items()), columns=['Number', 'pathway'])
df_legend=pd.merge(df_legend, unique_pathways_df_sorted[["pathway", "subSys"]], on='pathway', how="left")

heatmap_data_matched_filtered_2=heatmap_data_matched_filtered.copy()
mean_kd = heatmap_data_matched_filtered_2.mean(axis=1)
mean_resp = heatmap_data_matched_filtered_2.mean(axis=0)
mean_df_dict = {'Mean_kd': mean_kd,
                'Mean_resp': mean_resp}
mean_df = pd.DataFrame(mean_df_dict)
mean_df["pathway"]=mean_df.index
print(mean_df)
df_legend=pd.merge(df_legend, mean_df, on='pathway', how="left")
df_legend.to_excel("04_heatmap_legend.xlsx", "heatmap legend")
print("Made Number Legend for heatmap")

# Check if data to plot exists
if not heatmap_data_clipped.empty:
    # Heatmap Plot with gridspec for precise colorbar-placement
    plt.rcParams['font.size'] = 35/80*len(common_kept_pathways)
    fig = plt.figure(figsize=(max(10, heatmap_data_clipped.shape[1]), max(10, heatmap_data_clipped.shape[0])))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.05]) # Heatmap takes majority, colorbar smaller portion

    ax_heatmap = plt.subplot(gs[0])
    cbar_ax = plt.subplot(gs[1])
    width=1.0
    mapcolor="PiYG_r"
    axfontsize=35/80*len(common_kept_pathways)*3
    
    # Capture the colorbar object in a variable called `cbar`
    heatmap_plot = sns.heatmap(heatmap_data_clipped, cmap=mapcolor, center=0, annot=False,
                             linewidths=width, linecolor='black', square=True, ax=ax_heatmap, cbar_ax=cbar_ax)
    
    # Get the colorbar object from the heatmap plot and set the label and font size
    cbar = heatmap_plot.collections[0].colorbar
    cbar.set_label('Responses', fontdict={'fontsize': axfontsize})
    cbar.ax.tick_params(labelsize=axfontsize/3*2)

    # Make list of numbers of y ticks
    y_tick_labels = [pathway_number_mapping[pathway] for pathway in common_kept_pathways]

    # Set y-Ticks with numbers
    ax_heatmap.set_yticks(np.arange(len(common_kept_pathways)) + 0.5) # +0.5 centers labels
    ax_heatmap.set_yticklabels(y_tick_labels, rotation=0)

    # Same for x axis
    x_tick_labels = [pathway_number_mapping[pathway] for pathway in common_kept_pathways]
    ax_heatmap.set_xticks(np.arange(len(common_kept_pathways)) + 0.5)
    ax_heatmap.set_xticklabels(x_tick_labels, rotation=90)

    ax_heatmap.tick_params(axis='y', rotation=0)
    ax_heatmap.tick_params(axis='x', rotation=90)
    
    ax_heatmap.set_xlabel('Target-Pathway', fontsize=axfontsize)
    ax_heatmap.set_ylabel('Response-Pathway', fontsize=axfontsize)

    # Colorbar height correction
    # Make sure that heatmap was drawn
    fig.canvas.draw()

    # Set borders of heatmap data in data coordinates
    y_min_data, y_max_data = 83/95*heatmap_data_clipped.shape[0], heatmap_data_clipped.shape[0]

    # Transform Data coordinates in Pixel coordinates of Heatmap-Axes
    y_min_pixels, y_max_pixels = ax_heatmap.transData.transform([(0, y_min_data), (0, y_max_data)])[:, 1]

    # Calculate colorbar height in pixels
    colorbar_height_pixels = y_max_pixels - y_min_pixels

    # Get current position of colorbar axis in figure coordinates (0 to 1)
    cbar_pos = cbar_ax.get_position()
    cbar_x, cbar_y, cbar_width, cbar_height_fig = cbar_pos.bounds

    # Transform Pixelheight in figure coordinates
    fig_height_pixels = fig.canvas.get_width_height()[1]
    colorbar_height_fig_new = colorbar_height_pixels / fig_height_pixels

    # Set new colorbar height, keep y-position and width
    cbar_ax.set_position([cbar_x, cbar_y, cbar_width, colorbar_height_fig_new])

    plt.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(os.path.join(output_dir, f"04_pathway_interactions_heatmap_{cutoff}_and_{co2}_{mini}_num.png"))
    plt.savefig(os.path.join(output_dir, f"04_pathway_interactions_heatmap_{cutoff}_and_{co2}_{mini}_num.svg"))
    plt.close(fig)

else:
    print(f"No relevant data after filtering found. No heatmap made.")

print(f"All results stored in {output_dir}")

ecocyc_gene_pathways = pathway_table.groupby('pathway')['geneName'].apply(list).to_dict()

# Load the comprehensive pathway interaction data
heatmap_legend_df = pd.read_excel("04_heatmap_legend.xlsx", "heatmap legend", index_col="pathway")

# Initialize a list to store the pathway interactions that meet the criteria
pathways_in_legend = []

# Iterate through each pathway (row) in the heatmap data
for response_pathway, row in heatmap_legend_df.iterrows():
    pathways_in_legend.append({
        "Response Pathway": response_pathway,
    })

# Convert the list of interactions into a pandas DataFrame
pathways_in_legend_df = pd.DataFrame(pathways_in_legend)

# Analyze pathway knockdowns and create scatterplots
pathway_legend_results = {}

# Iterate through pathways
for index, interaction in pathways_in_legend_df.iterrows():
    response_pathway_html = interaction['Response Pathway']

    # Use a combined name for the pathway
    response_pathway_name_decoded_for_dict = response_pathway_html
    # Filter the merged table for the specific pathway interaction
    pathway_df = merged_table[merged_table['pathway']==response_pathway_html]
    pathway_gene_names = set(pathway_df['geneName'].tolist())

    all_pathway_kds = set(pathway_df.columns.tolist())
    ecocyc_kd_genes_in_pathway = set(ecocyc_gene_pathways[response_pathway_html])
    pathway_kds = all_pathway_kds.intersection(ecocyc_kd_genes_in_pathway)
    num_response_genes = len(pathway_gene_names.intersection(fc_table['geneName']))
    num_kds = len(pathway_kds.intersection(fc_table.columns))

    coverage_resp = 0
    coverage_kd = 0

    # Resp pathway coverage
    ecocyc_resp_genes_in_pathway = set(ecocyc_gene_pathways[response_pathway_html])
    overlap_resp = pathway_gene_names.intersection(ecocyc_resp_genes_in_pathway)
    coverage_resp = len(overlap_resp) / len(ecocyc_resp_genes_in_pathway) if len(ecocyc_resp_genes_in_pathway) > 0 else 0

    # Kd pathway coverage
    overlap_kd = pathway_kds.intersection(ecocyc_kd_genes_in_pathway)
    coverage_kd = len(overlap_kd) / len(ecocyc_kd_genes_in_pathway) if len(ecocyc_kd_genes_in_pathway) > 0 else 0

    pathway_legend_results[response_pathway_name_decoded_for_dict] = []
    pathway_legend_results[response_pathway_name_decoded_for_dict].append({'resp_coverage': coverage_resp, 'kd_coverage': coverage_kd, 'resp_num': num_response_genes, 'kd_num': num_kds})

# Save the pathway analysis results to an Excel file
coverage_pathway_data = []
for pathway, data in pathway_legend_results.items():
    for item in data:
        item['pathway'] = pathway
        coverage_pathway_data.append(item)

pathway_legend_results_df_01 = pd.DataFrame(coverage_pathway_data)

merged_legend_df= pd.merge(df_legend, pathway_legend_results_df_01, on='pathway', how='outer')
merged_legend_df=merged_legend_df.sort_values(by='Number')
merged_legend_df2=merged_legend_df.reset_index(drop=True)
merged_legend_df2.to_excel("04_heatmap_legend.xlsx", "heatmap legend")