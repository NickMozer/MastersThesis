# Scatterplot of log2 foldchanges of response metabolites versus proteins of a pathway responding to knockdowns from a strongly interacting other pathway

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import html

def decode_html_pathway_name(pathway_name_html):
    """Decodes HTML entities and prepares text for plot titles and filenames."""
    decoded_name = html.unescape(pathway_name_html)
    plot_title = decoded_name

    # Handle superscript for plot title and filename
    plot_title_sup = ""
    filename_sup = ""
    parts = decoded_name.split('<sup>')
    plot_title_sup += parts[0]
    filename_sup += parts[0]
    for part in parts[1:]:
        sub_parts = part.split('</sup>', 1)
        if len(sub_parts) == 2:
            plot_title_sup += r'$^{' + sub_parts[0] + '}$' + sub_parts[1]
            filename_sup += '^' + sub_parts[0] + sub_parts[1]
        else:
            plot_title_sup += '<sup>' + sub_parts[0]
            filename_sup += '<sup>' + sub_parts[0]
    plot_title = plot_title_sup

    # Handle subscript for plot title and filename
    plot_title_sub = ""
    filename_sub = ""
    parts = plot_title.split('<sub>')
    plot_title_sub += parts[0]
    filename_sub += parts[0]
    for part in parts[1:]:
        sub_parts = part.split('</sub>', 1)
        if len(sub_parts) == 2:
            plot_title_sub += r'$_{' + sub_parts[0] + '}$' + sub_parts[1]
            filename_sub += '_' + sub_parts[0] + sub_parts[1]
        else:
            plot_title_sub += '<sub>' + sub_parts[0]
            filename_sub += '<sub>' + sub_parts[0]
    plot_title = plot_title_sub
    return plot_title.strip()

def decode_html_pathway_name_2(pathway_name_html):
    """Decodes HTML entities and prepares text with style markers for plot titles."""
    decoded_name = html.unescape(pathway_name_html)
    styled_parts = []
    current_style = 'normal'
    current_text = ''

    def add_styled_part():
        nonlocal current_text, styled_parts, current_style
        if current_text:
            styled_parts.append((current_text, current_style))
            current_text = ''

    i = 0
    while i < len(decoded_name):
        if decoded_name[i:i+3] == '<i>' or decoded_name[i:i+3] == '<I>':
            add_styled_part()
            current_style = 'italic'
            i += 3
        elif decoded_name[i:i+4] == '</i>' or decoded_name[i:i+4] == '</I>':
            add_styled_part()
            current_style = 'normal'
            i += 4
        else:
            current_text += decoded_name[i]
            i += 1

    add_styled_part()
    return styled_parts

def create_safe_filename(pathway_name_html):
    """Creates a safe filename from the HTML pathway name."""
    decoded_name = html.unescape(pathway_name_html)
    safe_name = decoded_name.replace('<i>', '').replace('</i>', '').replace('pathway', '').replace('Pathway', '')
    safe_name = safe_name.replace('<I>', '').replace('</I>', '')
    safe_name = safe_name.replace('<sup>', '^').replace('</sup>', '')
    safe_name = safe_name.replace('<sub>', '_').replace('</sub>', '')
    safe_name = safe_name.replace(' ', '_').replace('(', '').replace(')', '').replace(',', '').replace("'", "").replace("<", "").replace(">", "").replace("/", "").replace('\\', "").replace(":", "")
    
    if len(safe_name) > 70:
        safe_name1 = f"{safe_name[:70]}"
    else:
        safe_name1 = safe_name
    return safe_name1

# Log2 fold change cutoff for response relevance
log2fc_cutoff = 1.0
# Minimum fraction of responding genes for plotting
min_response_fraction = 0.10

# Initialize date and output path
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_excel_path = "./sorted_by_path/26_pathway_fcs.xlsx"
output_path = os.path.join(output_dir, output_excel_path)

output_dir_very_relevant = "./sorted_by_path/very_relevant_0.3"
output_dir_relevant = "./sorted_by_path/relevant_0.3"
output_dir_less_relevant = "./sorted_by_path/less_relevant_0.3"
output_dir_not_relevant = "./sorted_by_path/not_relevant_0.3"
os.makedirs(output_dir_very_relevant, exist_ok=True)
os.makedirs(output_dir_relevant, exist_ok=True)
os.makedirs(output_dir_less_relevant, exist_ok=True)
os.makedirs(output_dir_not_relevant, exist_ok=True)

# Load fc, pathway tables
fc_table_prot = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
fc_table_metab = pd.read_excel("./219_metabolome_fc_table_low_rsds.xlsx", "Sheet1")
pathway_table = pd.read_excel("./20250327_pathways_reactants_ecoli.xlsx", "subpathways")

# Rename columns in pathway_table for metabolites
pathway_table = pathway_table.rename(columns={'metAbb': 'metabolite', 'Pathway': 'pathway', 'GeneName': 'geneName'})
fc_table_metab = fc_table_metab.rename(columns={'abbr': 'metabolite'})
# Get metabolite pathways for coverage calculation
ecocyc_metab_pathways = pathway_table.groupby('pathway')['metabolite'].apply(list).to_dict()
ecocyc_prot_pathways = pathway_table.groupby('pathway')['geneName'].apply(list).to_dict()
# Merge fc and pathway tables for metabolites
merged_table_metab = pd.merge(fc_table_metab, pathway_table[['metabolite', 'pathway']], on='metabolite', how='inner')
merged_table_prot = pd.merge(fc_table_prot, pathway_table[['geneName', 'pathway']], on='geneName', how='inner')

# Get the list of knockdown genes (columns in fc_table excluding 'metabolite')
knockdown_genes = [col for col in fc_table_metab.columns if col != 'metabolite']

# Load the list of pathway interactions with a response fraction of 0.3
interactions_03_df = pd.read_excel("26_pathway_interactions_at_0.3.xlsx", "Interactions_at_0.3")

# List to safe data of relevant plots
relevant_plot_data = []

# Analyze pathway knockdowns and create scatterplots
pathway_results_metab = {}
pathway_results_prot = {}
plotted_pathways_relevant = set()
plotted_pathways_not_relevant = set()
pathway_response_fractions_metab = {}
pathway_response_fractions_prot = {}

# Iterate through pathway interactions with a response fraction of 0.3
for index, interaction in interactions_03_df.iterrows():
    response_pathway_html = interaction['Response Pathway']
    knockdown_pathway_html = interaction['Knockdown Pathway']

    # Skip if it is the same pathway
    if response_pathway_html==knockdown_pathway_html:
        continue
    # Use a combined name for the pathway
    combined_pathway_name_html = f"{response_pathway_html} vs {knockdown_pathway_html}"
    styled_title_parts = decode_html_pathway_name_2(combined_pathway_name_html)
    safe_filename = create_safe_filename(combined_pathway_name_html)
    combined_pathway_name_decoded_for_dict = html.unescape(combined_pathway_name_html).replace('<i>', '').replace('</i>', '').replace('<I>', '').replace('</I>', '').replace('<sup>', '^').replace('</sup>', '').replace('<sub>', '_').replace('</sub>', '')

    # Filter the merged table for the specific pathway interaction
    pathway_df_metab = merged_table_metab[merged_table_metab['pathway']==response_pathway_html]
    pathway_metab_names = set(pathway_df_metab['metabolite'].tolist())

    all_pathway_kds_metab = set(pathway_df_metab.columns.tolist())
    ecocyc_kd_genes_in_pathway = set(ecocyc_prot_pathways[knockdown_pathway_html])
    pathway_kds = all_pathway_kds_metab.intersection(ecocyc_kd_genes_in_pathway)
    num_response_metabs = len(pathway_metab_names.intersection(fc_table_metab['metabolite']))
    num_kds = len(pathway_kds.intersection(fc_table_metab.columns))

    # Filter the merged table for the specific pathway interaction
    pathway_df_prot = merged_table_prot[merged_table_prot['pathway']==response_pathway_html]
    pathway_prot_names = set(pathway_df_prot['geneName'].tolist())

    num_response_prots = len(pathway_prot_names.intersection(fc_table_prot['geneName']))
    num_kds = len(pathway_kds.intersection(fc_table_prot.columns))

    # Condition 1: Minimum 2 response proteins and target genes in pathway
    if num_response_metabs < 2 or num_kds < 2 or num_response_prots < 2:
        continue

    # Find out coverage
    coverage_resp_metab = 0
    coverage_resp_prot = 0
    coverage_kd = 0

    # Resp pathway coverage metabs
    ecocyc_resp_metab_in_pathway = set(ecocyc_metab_pathways[response_pathway_html])
    overlap_resp_metab = pathway_metab_names.intersection(ecocyc_resp_metab_in_pathway)
    coverage_resp_metab = len(overlap_resp_metab) / len(ecocyc_resp_metab_in_pathway) if len(ecocyc_resp_metab_in_pathway) > 0 else 0

    # Resp pathway coverage prots
    ecocyc_resp_prot_in_pathway = set(ecocyc_prot_pathways[response_pathway_html])
    overlap_resp_prot = pathway_prot_names.intersection(ecocyc_resp_prot_in_pathway)
    coverage_resp_prot = len(overlap_resp_prot) / len(ecocyc_resp_prot_in_pathway) if len(ecocyc_resp_prot_in_pathway) > 0 else 0

    # Kd pathway coverage
    overlap_kd = pathway_kds.intersection(ecocyc_kd_genes_in_pathway)
    coverage_kd = len(overlap_kd) / len(ecocyc_kd_genes_in_pathway) if len(ecocyc_kd_genes_in_pathway) > 0 else 0

    if coverage_kd < 0.2 or coverage_resp_metab < 0.2 or coverage_resp_prot < 0.2:
        continue

    # Calculate response fraction for metabs
    pathway_results_metab[combined_pathway_name_decoded_for_dict] = []
    relevant_responses_count = 0
    total_responses = 0

    for response_metab in pathway_metab_names:
        if response_metab in fc_table_metab['metabolite'].values:
            for knocked_down_gene in knockdown_genes:
                if knocked_down_gene in pathway_kds:
                    try:
                        response_fc = fc_table_metab[fc_table_metab['metabolite'] == response_metab][knocked_down_gene].values[0]
                        pathway_results_metab[combined_pathway_name_decoded_for_dict].append({'response_metab': response_metab, 'knockdown_gene': knocked_down_gene,  'response_fc': response_fc, 'in_pathway': True})
                        if response_fc >= log2fc_cutoff and knocked_down_gene != response_metab:
                            relevant_responses_count += 1
                        if knocked_down_gene != response_metab:
                            total_responses += 1
                    except IndexError:
                        print(f"response metab '{response_metab}' not found in fc_table_metab for knockdown of '{knocked_down_gene}' in pathway '{combined_pathway_name_html}'.")
                else:
                    try:
                        other_response_fc = fc_table_metab[fc_table_metab['metabolite'] == response_metab][knocked_down_gene].values[0]
                        pathway_results_metab[combined_pathway_name_decoded_for_dict].append({'response_metab': response_metab, 'knockdown_gene': knocked_down_gene,  'response_fc': other_response_fc, 'in_pathway': False})
                    except IndexError:
                        pass # Gene might not be present for this knockdown

    response_fraction_metab = relevant_responses_count / total_responses if total_responses > 0 else 0
    pathway_response_fractions_metab[combined_pathway_name_decoded_for_dict] = response_fraction_metab

    # Calculate response fraction for prots
    pathway_results_prot[combined_pathway_name_decoded_for_dict] = []
    relevant_responses_count = 0
    total_responses = 0

    for response_prot in pathway_prot_names:
        if response_prot in fc_table_prot['geneName'].values:
            for knocked_down_gene in knockdown_genes:
                if knocked_down_gene in pathway_kds:
                    try:
                        response_fc = fc_table_prot[fc_table_prot['geneName'] == response_prot][knocked_down_gene].values[0]
                        pathway_results_prot[combined_pathway_name_decoded_for_dict].append({'response_prot': response_prot, 'knockdown_gene': knocked_down_gene,  'response_fc': response_fc, 'in_pathway': True})
                        if response_fc >= log2fc_cutoff and knocked_down_gene != response_prot:
                            relevant_responses_count += 1
                        if knocked_down_gene != response_prot:
                            total_responses += 1
                    except IndexError:
                        print(f"response prot '{response_prot}' not found in fc_table_prot for knockdown of '{knocked_down_gene}' in pathway '{combined_pathway_name_html}'.")
                else:
                    try:
                        other_response_fc = fc_table_prot[fc_table_prot['geneName'] == response_prot][knocked_down_gene].values[0]
                        pathway_results_prot[combined_pathway_name_decoded_for_dict].append({'response_prot': response_prot, 'knockdown_gene': knocked_down_gene,  'response_fc': other_response_fc, 'in_pathway': False})
                    except IndexError:
                        pass # Gene might not be present for this knockdown

    response_fraction_prot = relevant_responses_count / total_responses if total_responses > 0 else 0
    pathway_response_fractions_prot[combined_pathway_name_decoded_for_dict] = response_fraction_prot

    # Create and save scatterplot
    if pathway_results_metab[combined_pathway_name_decoded_for_dict] and pathway_results_prot[combined_pathway_name_decoded_for_dict]:
        results_df_metab = pd.DataFrame(pathway_results_metab[combined_pathway_name_decoded_for_dict])
        results_df_prot = pd.DataFrame(pathway_results_prot[combined_pathway_name_decoded_for_dict])
        plt.rcParams['font.size'] = 18
        fig = plt.figure(figsize=(8, 10))

        # sort data for other response genes
        pathway_metab_responses_df = results_df_metab[results_df_metab['in_pathway'] == True].copy()
        pathway_prot_responses_df = results_df_prot[results_df_prot['in_pathway'] == True].copy()

        y_max_metab = pathway_metab_responses_df['response_fc'].max()
        y_min_metab = pathway_metab_responses_df['response_fc'].min()
        y_max_prot = pathway_prot_responses_df['response_fc'].max()
        y_min_prot = pathway_prot_responses_df['response_fc'].min()
        global_max = max(y_max_metab, y_max_prot)
        global_min = min(y_min_metab, y_min_prot)
        global_max=global_max*1.05
        global_min=global_min*1.05

        pathway_kds_list=list(pathway_kds)
        kd_color_mapping={}

        colors=("teal", "blue", "firebrick", "deeppink", "navy", "peru", "purple", "grey", "black", "green")
        colors2=("navy", "peru", "purple", "darkgrey", "black", "green")
        # specifically color pentose pathways differently 
        if "pentose" in safe_filename: 
            colors=colors2
        for i, current_color in enumerate(colors):
            if num_kds > i:
                kd_color_mapping[pathway_kds_list[i]]=current_color
        
        # Add 'kd_color' column to pathway_metab_responses_df
        pathway_metab_responses_df['kd_color'] = pathway_metab_responses_df['knockdown_gene'].map(kd_color_mapping)

        # Add 'kd_color' column to pathway_prot_responses_df
        pathway_prot_responses_df['kd_color'] = pathway_prot_responses_df['knockdown_gene'].map(kd_color_mapping)

        # Scatter Plot for response metabolites (left side)
        texts_prot = []
        labeled_points_pathway = set()
        texts_metab = []
        labeled_points_pathway = set()
        s1=7.6
        ax1 = plt.subplot(1, 2, 1) # 1 row, 2 columns, first Subplot
        sns.stripplot(y="response_fc", data=pathway_metab_responses_df, hue="knockdown_gene", size=s1, jitter=False, ax=ax1, palette=kd_color_mapping)
        ax1.set_ylim(global_min, global_max)
        ax1.set_xticks([])
        ax1.set_ylabel('Log2 fold change')
        ax1.set_xlabel('Response\nmetabolites')
        ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        for pos in ['right', 'top']:
            ax1.spines[pos].set_visible(False)

        # Adjust position of metab labels
        for resp_gene in pathway_metab_responses_df['response_metab'].unique():
            subset_pathway = pathway_metab_responses_df[(pathway_metab_responses_df['response_metab'] == resp_gene) & (pathway_metab_responses_df['in_pathway'] == True)].sort_values('response_fc')
            n_points_pathway = len(subset_pathway)
            labels_to_add_pathway = []

            for i in range(n_points_pathway):
                labels_to_add_pathway.append(subset_pathway.iloc[i])
            for row in labels_to_add_pathway:
                point_tuple = (row['response_metab'], row['knockdown_gene'])
                if point_tuple not in labeled_points_pathway:
                    text_metab=ax1.text(x=0.002, y=row['response_fc'], s=row['response_metab'], color=row['kd_color'])
                    texts_metab.append(text_metab)
                    labeled_points_pathway.add(point_tuple)
        adjust_text(texts_metab, arrowprops=None, autoalign='xy', only_move={'points':'x', 'text':'x'})

        ax2 = plt.subplot(1, 2, 2, sharey=ax1) # 1 row, 2 columns, second Subplot, shares y-axis
        sns.stripplot(y="response_fc", data=pathway_prot_responses_df, hue="knockdown_gene", size=s1, jitter=False, ax=ax2, palette=kd_color_mapping)
        ax2.set_xticks([])
        ax2.set_xlabel('Response\nproteins')
        ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
        ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax2.set_ylabel('') # No new y-label
        for pos in ['left', 'top', 'right']:
            ax2.spines[pos].set_visible(False)

        # Adjust position of labels
        for resp_gene in pathway_prot_responses_df['response_prot'].unique():
            subset_pathway = pathway_prot_responses_df[(pathway_prot_responses_df['response_prot'] == resp_gene) & (pathway_prot_responses_df['in_pathway'] == True)].sort_values('response_fc')
            n_points_pathway = len(subset_pathway)
            labels_to_add_pathway = []

            for i in range(n_points_pathway):
                labels_to_add_pathway.append(subset_pathway.iloc[i])
            for row in labels_to_add_pathway:
                point_tuple = (row['response_prot'], row['knockdown_gene'])
                if point_tuple not in labeled_points_pathway:
                    text_prot=ax2.text(x=0.002, y=row['response_fc'], s=row['response_prot'], color=row['kd_color'])
                    texts_prot.append(text_prot)
                    labeled_points_pathway.add(point_tuple)

        adjust_text(texts_prot, arrowprops=None, autoalign='xy', only_move={'points':'x', 'text':'x'})

        leg1=ax1.legend()
        leg1.set_visible(False) # Make this legend invisible
        leg2=ax2.legend(title='Target', bbox_to_anchor=(0.95, 1), loc='upper left')
        leg2.set_visible(True) # Make this legend visible

        # Format title
        title_parts = []
        for text, style in styled_title_parts:
            if style == 'italic':
                title_parts.append(f'$\\it{{{text}}}$')
            else:
                title_parts.append(text)
        title_with_italics=" ".join(title_parts)
        title_with_sub_and_sup=decode_html_pathway_name(title_with_italics)
        title_with_sub_and_sup = title_with_sub_and_sup.replace('  ', ' ').replace('( ', '(').replace(' -', '-').replace(' )', ')').replace(f'$\\it{{denovo}}$', f'$\\it{{de novo}}$').replace(f'$\\it{{E.coli}}$', f'$\\it{{E. coli}}$')
        woerter = title_with_sub_and_sup.split()

        if woerter:
            erstes_wort = woerter[0]
            kapitalisiertes_erstes_wort = erstes_wort.capitalize()
            woerter[0] = kapitalisiertes_erstes_wort
            title_with_sub_and_sup = " ".join(woerter)

        formatted_title = f'{title_with_sub_and_sup}\n(Response metabolites: {response_fraction_metab:.2f}, Response proteins: {response_fraction_prot:.2f},\nCoverage_kd: {coverage_kd:.2f}), Coverage_resp_metab: {coverage_resp_metab:.2f}, Coverage_resp_prot: {coverage_resp_prot:.2f})'

        fig.tight_layout()
        plt.rcParams['svg.fonttype'] = 'none'

        # Save in two folders to prefilter more relevant pathways
        if response_fraction_metab >= 0.5 and num_kds >= 3 and num_response_metabs >= 3 and coverage_kd >= 0.90 and coverage_resp_metab >= 0.90:
            plt.savefig(os.path.join(output_dir_very_relevant, f"26_{safe_filename}_colored3.svg"))
            plt.savefig(os.path.join(output_dir_very_relevant, f"26_{safe_filename}_colored3.png"))
            plotted_pathways_relevant.add(combined_pathway_name_decoded_for_dict)
        elif response_fraction_metab >= 0.3 and num_kds >= 3 and num_response_metabs >= 3 and coverage_kd >= 0.70 and coverage_resp_metab >= 0.70:
            plt.savefig(os.path.join(output_dir_relevant, f"26_{safe_filename}_colored3.svg"))
            plt.savefig(os.path.join(output_dir_relevant, f"26_{safe_filename}_colored3.png"))
            plotted_pathways_relevant.add(combined_pathway_name_decoded_for_dict)
        elif response_fraction_metab >= 0.1 and num_kds >= 3 and num_response_metabs >= 3 and coverage_kd >= 0.50 and coverage_resp_metab >= 0.50:
            plt.savefig(os.path.join(output_dir_less_relevant, f"26_{safe_filename}_colored3.svg"))
            plt.savefig(os.path.join(output_dir_less_relevant, f"26_{safe_filename}_colored3.png"))
            plotted_pathways_relevant.add(combined_pathway_name_decoded_for_dict)
        else:
            plt.savefig(os.path.join(output_dir_not_relevant, f"26_{safe_filename}_colored3.svg"))
            plt.savefig(os.path.join(output_dir_not_relevant, f"26_{safe_filename}_colored3.png"))
            plotted_pathways_not_relevant.add(combined_pathway_name_decoded_for_dict)
        plt.close()

print(f"Scatterplots for {len(plotted_pathways_relevant)} pathways (meeting the criteria) saved to: {output_dir_relevant}")
print(f"Scatterplots for {len(plotted_pathways_not_relevant)} pathways (not meeting the criteria) saved to: {output_dir_not_relevant}")

# Save the pathway analysis results to an Excel file
relevant_pathway_data_metab = []
irrelevant_pathway_data_metab = []
for pathway, data in pathway_results_metab.items():
    if pathway in plotted_pathways_relevant:
        for item in data:
            item['pathway'] = pathway
            relevant_pathway_data_metab.append(item)
    elif pathway in plotted_pathways_not_relevant:
        for item in data:
            item['pathway'] = pathway
            irrelevant_pathway_data_metab.append(item)

if relevant_pathway_data_metab or irrelevant_pathway_data_metab:
    pathway_fraction_df = pd.DataFrame(list(pathway_response_fractions_metab.items()), columns=['pathway', f'response_fraction_log2fc_{log2fc_cutoff}'])
    relevant_pathways_df = pd.DataFrame(sorted(list(plotted_pathways_relevant)), columns=['relevant_pathway_name'])
    if relevant_pathway_data_metab:
        pathway_results_df_01 = pd.DataFrame(relevant_pathway_data_metab)
        merged_output_df_relevant = pd.merge(pathway_results_df_01, pathway_fraction_df, on='pathway', how='left')
        merged_output_df_relevant_densed=merged_output_df_relevant[["pathway", "response_fraction_log2fc_1.0"]].drop_duplicates()
    if irrelevant_pathway_data_metab:
        pathway_results_df_02 = pd.DataFrame(irrelevant_pathway_data_metab)
        merged_output_df_irrelevant = pd.merge(pathway_results_df_02, pathway_fraction_df, on='pathway', how='left')
        merged_output_df_irrelevant_densed=merged_output_df_irrelevant[["pathway", "response_fraction_log2fc_1.0"]].drop_duplicates()
    with pd.ExcelWriter(output_path) as writer:
        if relevant_pathway_data_metab:
            merged_output_df_relevant_densed.to_excel(writer, sheet_name="relevant_pathway_metab_fcs", index=False)
        if irrelevant_pathway_data_metab:
            merged_output_df_irrelevant_densed.to_excel(writer, sheet_name="irrelevant_pathway_metab_fcs", index=False)
        relevant_pathways_df.to_excel(writer, sheet_name="relevant_pathway_names", index=False)
    print(f"Pathway analysis results for all analyzed pathways saved to: {output_path}")
else:
    print("No pathway analysis results to save (no pathways met the initial criteria).")

# Save the pathway analysis results to an Excel file
relevant_pathway_data_prot = []
irrelevant_pathway_data_prot = []
for pathway, data in pathway_results_prot.items():
    if pathway in plotted_pathways_relevant:
        for item in data:
            item['pathway'] = pathway
            relevant_pathway_data_prot.append(item)
    elif pathway in plotted_pathways_not_relevant:
        for item in data:
            item['pathway'] = pathway
            irrelevant_pathway_data_prot.append(item)

if relevant_pathway_data_prot or irrelevant_pathway_data_prot:
    pathway_fraction_df = pd.DataFrame(list(pathway_response_fractions_prot.items()), columns=['pathway', f'response_fraction_log2fc_{log2fc_cutoff}'])
    relevant_pathways_df = pd.DataFrame(sorted(list(plotted_pathways_relevant)), columns=['relevant_pathway_name'])
    if relevant_pathway_data_prot:
        pathway_results_df_01 = pd.DataFrame(relevant_pathway_data_prot)
        merged_output_df_relevant = pd.merge(pathway_results_df_01, pathway_fraction_df, on='pathway', how='left')
        merged_output_df_relevant_densed=merged_output_df_relevant[["pathway", "response_fraction_log2fc_1.0"]].drop_duplicates()
    if irrelevant_pathway_data_prot:
        pathway_results_df_02 = pd.DataFrame(irrelevant_pathway_data_prot)
        merged_output_df_irrelevant = pd.merge(pathway_results_df_02, pathway_fraction_df, on='pathway', how='left')
        merged_output_df_irrelevant_densed=merged_output_df_irrelevant[["pathway", "response_fraction_log2fc_1.0"]].drop_duplicates()
    with pd.ExcelWriter(output_path) as writer:
        if relevant_pathway_data_prot:
            merged_output_df_relevant_densed.to_excel(writer, sheet_name="relevant_pathway_prot_fcs", index=False)
        if irrelevant_pathway_data_prot:
            merged_output_df_irrelevant_densed.to_excel(writer, sheet_name="irrelevant_pathway_prot_fcs", index=False)
        relevant_pathways_df.to_excel(writer, sheet_name="relevant_pathway_names", index=False)
    print(f"Pathway analysis results for all analyzed pathways saved to: {output_path}")
else:
    print("No pathway analysis results to save (no pathways met the initial criteria).")