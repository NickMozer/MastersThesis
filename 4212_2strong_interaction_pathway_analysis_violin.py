# Scatterplots of log2 foldchanges of response genes in a strongly interacting pathway of the target/knockdown gene (rsds below 50%)
# No legends, plot sizes adapted to number of knockdown/target genes
# Minimum 2 response genes per knockdown, minimum 20% coverage of the pathway
# 95% Coverage and minimum 3 resp and kd genes pathways are sorted into an extra folder
# Pathways that didnt meet that last condition get an extra folder of less relevance for themselves
# Decode html format of pathway names
# Add violins of other response genes in gray

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
output_excel_path = "./sorted_by_path/13_pathway_fcs.xlsx"
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
fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
pathway_table = pd.read_excel("./20250327_pathways_reactants_ecoli.xlsx", "subpathways")
pathway_table = pathway_table.rename(columns={'GeneName': 'geneName', 'Pathway': 'pathway'})
ecocyc_gene_pathways = pathway_table.groupby('pathway')['geneName'].apply(list).to_dict()

# Merge fc and pathway tables
merged_table = pd.merge(fc_table, pathway_table[['geneName', 'pathway']], on='geneName', how='inner')

# Get the list of knockdown genes (columns in fc_table excluding 'geneName')
knockdown_genes = [col for col in fc_table.columns if col != 'geneName']

# Load the list of pathway interactions with a response fraction of 0.3
interactions_03_df = pd.read_excel("13_pathway_interactions_at_0.3.xlsx", "Interactions_at_0.3")

# List to safe data of relevant plots
relevant_plot_data = []

# Analyze pathway knockdowns and create scatterplots
pathway_results = {}
plotted_pathways_relevant = set()
plotted_pathways_not_relevant = set()
pathway_response_fractions = {}

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
    pathway_df = merged_table[merged_table['pathway']==response_pathway_html]
    pathway_gene_names = set(pathway_df['geneName'].tolist())

    all_pathway_kds = set(pathway_df.columns.tolist())
    ecocyc_kd_genes_in_pathway = set(ecocyc_gene_pathways[knockdown_pathway_html])
    pathway_kds = all_pathway_kds.intersection(ecocyc_kd_genes_in_pathway)
    num_response_genes = len(pathway_gene_names.intersection(fc_table['geneName']))
    num_kds = len(pathway_kds.intersection(fc_table.columns))

    # Condition 1: Minimum 2 response genes in pathway
    if num_response_genes < 2 or num_kds < 2:
        continue

    # Condition 2: Minimum 20% coverage of Ecocyc-genes in pathway
    coverage_resp = 0
    coverage_kd = 0

    # Resp pathway coverage
    ecocyc_resp_genes_in_pathway = set(ecocyc_gene_pathways[response_pathway_html])
    overlap_resp = pathway_gene_names.intersection(ecocyc_resp_genes_in_pathway)
    coverage_resp = len(overlap_resp) / len(ecocyc_resp_genes_in_pathway) if len(ecocyc_resp_genes_in_pathway) > 0 else 0

    # Kd pathway coverage
    overlap_kd = pathway_kds.intersection(ecocyc_kd_genes_in_pathway)
    coverage_kd = len(overlap_kd) / len(ecocyc_kd_genes_in_pathway) if len(ecocyc_kd_genes_in_pathway) > 0 else 0

    if coverage_kd < 0.2 or coverage_resp < 0.2:
        continue

    # Calculate response fraction
    pathway_results[combined_pathway_name_decoded_for_dict] = []
    relevant_responses_count = 0
    total_responses = 0

    for response_gene in pathway_gene_names:
        if response_gene in fc_table['geneName'].values:
            for knocked_down_gene in knockdown_genes:
                if knocked_down_gene in pathway_kds:
                    try:
                        response_fc = fc_table[fc_table['geneName'] == response_gene][knocked_down_gene].values[0]
                        pathway_results[combined_pathway_name_decoded_for_dict].append({'response_gene': response_gene, 'knockdown_gene': knocked_down_gene,  'response_fc': response_fc, 'in_pathway': True})
                        if response_fc >= log2fc_cutoff and knocked_down_gene != response_gene:
                            relevant_responses_count += 1
                        if knocked_down_gene != response_gene:
                            total_responses += 1
                    except IndexError:
                        print(f"response gene '{response_gene}' not found in fc_table for knockdown of '{knocked_down_gene}' in pathway '{combined_pathway_name_html}'.")
                else:
                    try:
                        other_response_fc = fc_table[fc_table['geneName'] == response_gene][knocked_down_gene].values[0]
                        pathway_results[combined_pathway_name_decoded_for_dict].append({'response_gene': response_gene, 'knockdown_gene': knocked_down_gene,  'response_fc': other_response_fc, 'in_pathway': False})
                    except IndexError:
                        pass # Gene might not be present for this knockdown

    response_fraction = relevant_responses_count / total_responses if total_responses > 0 else 0
    pathway_response_fractions[combined_pathway_name_decoded_for_dict] = response_fraction

    # Create and save scatterplot
    if pathway_results[combined_pathway_name_decoded_for_dict]:
        results_df = pd.DataFrame(pathway_results[combined_pathway_name_decoded_for_dict])

        # Extract order of response genes from pathway table
        current_pathway_genes = pathway_table[pathway_table['pathway']==response_pathway_html]['geneName'].unique().tolist()

        # Create mapping for sorting order 
        gene_order_mapping = pd.Series(range(len(current_pathway_genes)), index=current_pathway_genes)

        # Use .categorical for sorting
        results_df['response_gene_categorical'] = pd.Categorical(results_df['response_gene'], categories=current_pathway_genes, ordered=True)

        # Sort results_df for categorial column
        results_df = results_df.sort_values(by='response_gene_categorical')

        # Discard extra column
        results_df = results_df.drop(columns=['response_gene_categorical'])

        response_gene_names = results_df['response_gene'].unique().tolist()
        # List of italics names
        italic_labels = [f'$\\it{{{gene}}}$' for gene in response_gene_names]

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1)

        # sort data for other response genes
        other_knockdowns_df = results_df[results_df['in_pathway'] == False]
        pathway_responses_df = results_df[results_df['in_pathway'] == True]
        significant_other_responses = other_knockdowns_df[abs(other_knockdowns_df['response_fc']) >= 1.5]
        non_significant_other_responses = other_knockdowns_df[abs(other_knockdowns_df['response_fc']) < 1.5]

        # Violin width relative to fogure width
        num_resp = len(response_gene_names)
        if num_resp == 1:
            violin_width = 0.1
            # Make explicit x-ticks and limits for enough space
            ax.set_xticks([0])
            ax.set_xticklabels(italic_labels)
            ax.set_xlim(-0.75, 0.75) # Space left and right
        elif len(response_gene_names) == 2:
            violin_width = 0.15
            positions = [0 - violin_width/2 - 0.1, 0 + violin_width/2 + 0.1] # Adjust positions
            ax.set_xticklabels(italic_labels)
            ax.set_xlim(positions[0] - 0.25, positions[1] + 0.25)
        else:
            if len(response_gene_names) == 3:
                violin_width = 0.2
            elif len(response_gene_names) == 4:
                violin_width = 0.25
            elif len(response_gene_names) == 5:
                violin_width = 0.3
            else:
                violin_width = 0.5
            # Set x-ticks based on number of genes
            ax.set_xticks(range(len(response_gene_names)))
            # Use italics list for ax labels
            ax.set_xticklabels(italic_labels)
        # Violin Plot for non significant other responses
        sns.violinplot(x='response_gene', y='response_fc', data=non_significant_other_responses,
                            color='darkgray', alpha=0.5, inner='quartile', linewidth=0.5, label='Other Responses', zorder=1, width=violin_width)

        # Scatterplot for significant other Response Genes
        plt.scatter(x=significant_other_responses['response_gene'], y=significant_other_responses['response_fc'],
                            color='darkgray', alpha=0.5, s=20, edgecolors='darkgray', label='Other Responses (|log2FC| >= 1.5)', facecolors='lightgray',zorder=2)

        # Plot pathway responses (dodgerblue)
        plt.scatter(x=pathway_responses_df['response_gene'], y=pathway_responses_df['response_fc'],
                            color='dodgerblue', alpha=0.5, linewidths=0.5, edgecolors='dodgerblue', label='Pathway Responses', facecolors='dodgerblue', zorder=3)

        plt.xlabel('Response gene')
        plt.ylabel(f'Log2 fold change')

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

        formatted_title = f'{title_with_sub_and_sup}\n(Response: {response_fraction:.2f}, Coverage_kd: {coverage_kd:.2f}), Coverage_resp: {coverage_resp:.2f})'
        plt.title(formatted_title)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
        plt.tick_params(axis='y', which='both', right=False, left=True, labelleft=True)
        for pos in ['right', 'top', 'bottom']:
            plt.gca().spines[pos].set_visible(False)
        plt.xticks(rotation=45, ha='right')
        leg=plt.legend(title='Response Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        leg.set_visible(False) # Make the legend invisible

        # Adjust position of labels
        texts = []
        labeled_points_pathway = set()
        for resp_gene in results_df['response_gene'].unique():
            subset_pathway = results_df[(results_df['response_gene'] == resp_gene) & (results_df['in_pathway'] == True)].sort_values('response_fc')
            n_points_pathway = len(subset_pathway)
            labels_to_add_pathway = []

            for i in range(n_points_pathway):
                labels_to_add_pathway.append(subset_pathway.iloc[i])
            for row in labels_to_add_pathway:
                point_tuple = (row['response_gene'], row['knockdown_gene'])
                if point_tuple not in labeled_points_pathway:
                    texts.append(plt.text(row['response_gene'], row['response_fc'], row['knockdown_gene'], color='dodgerblue'))
                    labeled_points_pathway.add(point_tuple)

        adjust_text(texts, arrowprops=None, autoalign='xy', only_move={'points':'xy', 'text':'xy'})

        plt.tight_layout()
        plt.rcParams['svg.fonttype'] = 'none'

        # Save in seperate folders to prefilter more relevant pathways
        if response_fraction >= 0.5 and num_kds >= 3 and num_response_genes >= 3 and coverage_kd >= 0.90 and coverage_resp >= 0.90:
            plt.savefig(os.path.join(output_dir_very_relevant, f"13_sorted_by_path_{safe_filename}.svg"))
            plt.savefig(os.path.join(output_dir_very_relevant, f"13_sorted_by_path_{safe_filename}.png"))
            plotted_pathways_relevant.add(combined_pathway_name_decoded_for_dict)
        elif response_fraction >= 0.3 and num_kds >= 3 and num_response_genes >= 3 and coverage_kd >= 0.70 and coverage_resp >= 0.70:
            plt.savefig(os.path.join(output_dir_relevant, f"13_sorted_by_path_{safe_filename}.svg"))
            plt.savefig(os.path.join(output_dir_relevant, f"13_sorted_by_path_{safe_filename}.png"))
            plotted_pathways_relevant.add(combined_pathway_name_decoded_for_dict)
        elif response_fraction >= 0.1 and num_kds >= 3 and num_response_genes >= 3 and coverage_kd >= 0.50 and coverage_resp >= 0.50:
            plt.savefig(os.path.join(output_dir_less_relevant, f"13_sorted_by_path_{safe_filename}.svg"))
            plt.savefig(os.path.join(output_dir_less_relevant, f"13_sorted_by_path_{safe_filename}.png"))
            plotted_pathways_relevant.add(combined_pathway_name_decoded_for_dict)
        else:
            plt.savefig(os.path.join(output_dir_not_relevant, f"13_sorted_by_path_{safe_filename}.svg"))
            plt.savefig(os.path.join(output_dir_not_relevant, f"13_sorted_by_path_{safe_filename}.png"))
            plotted_pathways_not_relevant.add(combined_pathway_name_decoded_for_dict)
        plt.close()

print(f"Scatterplots for {len(plotted_pathways_relevant)} pathways (meeting the criteria) saved to: {output_dir_relevant}")
print(f"Scatterplots for {len(plotted_pathways_not_relevant)} pathways (not meeting the criteria) saved to: {output_dir_not_relevant}")

# Save the pathway analysis results to an Excel file
relevant_pathway_data = []
irrelevant_pathway_data = []
for pathway, data in pathway_results.items():
    if pathway in plotted_pathways_relevant:
        for item in data:
            item['pathway'] = pathway
            relevant_pathway_data.append(item)
    elif pathway in plotted_pathways_not_relevant:
        for item in data:
            item['pathway'] = pathway
            irrelevant_pathway_data.append(item)

if relevant_pathway_data or irrelevant_pathway_data:
    pathway_fraction_df = pd.DataFrame(list(pathway_response_fractions.items()), columns=['pathway', f'response_fraction_log2fc_{log2fc_cutoff}'])
    relevant_pathways_df = pd.DataFrame(sorted(list(plotted_pathways_relevant)), columns=['relevant_pathway_name'])
    if relevant_pathway_data:
        pathway_results_df_01 = pd.DataFrame(relevant_pathway_data)
        merged_output_df_relevant = pd.merge(pathway_results_df_01, pathway_fraction_df, on='pathway', how='left')
        merged_output_df_relevant_densed=merged_output_df_relevant[["pathway", "response_fraction_log2fc_1.0"]].drop_duplicates()
    if irrelevant_pathway_data:
        pathway_results_df_02 = pd.DataFrame(irrelevant_pathway_data)
        merged_output_df_irrelevant = pd.merge(pathway_results_df_02, pathway_fraction_df, on='pathway', how='left')
        merged_output_df_irrelevant_densed=merged_output_df_irrelevant[["pathway", "response_fraction_log2fc_1.0"]].drop_duplicates()
    with pd.ExcelWriter(output_path) as writer:
        if relevant_pathway_data:
            merged_output_df_relevant_densed.to_excel(writer, sheet_name="relevant_pathway_knockdown_fcs", index=False)
        if irrelevant_pathway_data:
            merged_output_df_irrelevant_densed.to_excel(writer, sheet_name="irrelevant_pathway_knockdown_fcs", index=False)
        relevant_pathways_df.to_excel(writer, sheet_name="relevant_pathway_names", index=False)
    print(f"Pathway analysis results for all analyzed pathways saved to: {output_path}")
else:
    print("No pathway analysis results to save (no pathways met the initial criteria).")
