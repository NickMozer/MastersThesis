# Scatter plot the CRISPR knockdown strains Ale proteome fold change data against the corresponding Donati strains fold change data
# Same axes within plot and throughout all plots (highest log2 fc were never above 8, so 8 was chosen as the border)
# using final filtered knockdown table
# black points, increased size to 512
# include pearson R and diagonal line
# Take all svgs made in this code and make one 4x5 svg figure from it (export png in Inkscape after)
# Removing of single x and y labels in plot for all single figures, add common labels

# Import libraries
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import xml.etree.ElementTree as ET

f_size=56
title_font_size = "56"
s1=512

cols = 5
rows = 4

common_x_title = "Log2 fold change dataset 1"
common_y_title = "Log2 fold change dataset 2"
axes_offset = 75 # Space for ax titles

# Initialize output path
output_excel_path = "01_dataset_comparison.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_dir0 = os.path.join(output_dir, './01_single_comparison_figures')
os.makedirs(output_dir0, exist_ok=True)
output_path = os.path.join(output_dir, output_excel_path)
output_filename = f"01_dataset_comparisons_h2_{rows}_{cols}.svg"

# Load donati, fc tables
donati_table = pd.read_excel("./1-s2.0-S240547122030418X-mmc9.xlsx", "Table 7")
fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")

# Remove first row (two general column heads)
donati_table.columns=donati_table.iloc[0]
donati_table=donati_table.iloc[1:]

# Pick Donati fcs
donati_table_fcs=donati_table.iloc[:, :31]
donati_table_fcs=donati_table_fcs[
    list(donati_table_fcs.columns[0:1]) + list(donati_table_fcs.columns[2:31])
]
donati_table_fcs = donati_table_fcs.rename(columns={'Gene Name': 'geneName'})

# merge fcs, donati fcs and donati errors on geneName, ensure number data format
for col in donati_table_fcs.columns[1:]: 
    donati_table_fcs[col] = pd.to_numeric(donati_table_fcs[col], errors='coerce')

for col in fc_table.columns[1:]: 
    fc_table[col] = pd.to_numeric(fc_table[col], errors='coerce')
comparison_table = pd.merge(donati_table_fcs, fc_table, on='geneName', how='inner', suffixes=('_01', '_02'))

print(comparison_table.head())
comparison_table=comparison_table.loc[:, ["geneName"] + [col for col in comparison_table.columns if "_0" in col]]
comparison_table = comparison_table[list(comparison_table.columns[:1]) + sorted(comparison_table.columns[1:])]
print(comparison_table.head())

column_names = []
pearson_r_values = []

# Check which Donati CRISPR strains are also in new dataset
for column in donati_table_fcs.columns:
    if "gene" not in column and column in fc_table.columns:
        
        # Scatter plot of fcs
        x=comparison_table[f"{column}_01"]
        y=comparison_table[f"{column}_02"]
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(x, y, 
                    s=s1,
                    marker='o',
                    edgecolor='black',
                    facecolor='black',
                    alpha=1.0)

        # Set global minimal- and maximum values for axes, as maximum values were never above or below 8 for the comparisons before
        global_min = -8
        global_max = 8

        # Get the linewidth of one of the axes spines
        axis_linewidth = ax.spines['bottom'].get_linewidth()
        axis_linewidth_4=axis_linewidth*4

        # Add diagonal line with the same linewidth as the axes
        ax.plot([global_min, global_max], [global_min, global_max], color='black', linestyle='-', linewidth=axis_linewidth_4)

        # Calculate Pearson correlation coefficient
        correlation, _ = pearsonr(x, y)
        ax.text(global_min + 0.5, global_max - 1.5, f'R: {correlation:.2f}', fontsize=f_size)

        column_names.append(column)
        pearson_r_values.append(correlation)

        # Determine ax ranges
        ax.set_xlim(global_min, global_max)
        ax.set_ylim(global_min, global_max)

        ax.set_title(f"$\\it{{{column}}}$", fontsize=f_size)
        ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
        ax.tick_params(axis='y', which='both', right=False, left=True, labelleft=True)
        for pos in ['right', 'top']:
            ax.spines[pos].set_visible(False)
        for pos in ['left', 'bottom']:
            ax.spines[pos].set_linewidth(axis_linewidth_4)

        # Set unique IDs for the x and y axes
        ax.xaxis.set_label_text("Donati Fold Change")
        ax.yaxis.set_label_text("New Dataset Fold Change")
        ax.xaxis.label.set_gid("x_axis_" + column)
        ax.yaxis.label.set_gid("y_axis_" + column)
        gewuenschte_ticks = [-8, -4, 0, 4, 8]
        ax.set_xticks(gewuenschte_ticks)
        ax.set_yticks(gewuenschte_ticks)
        ax.xaxis.set_ticklabels(ax.get_xticklabels(), gid="x_axis_ticks_" + column, fontsize=f_size)
        ax.yaxis.set_ticklabels(ax.get_yticklabels(), gid="y_axis_ticks_" + column, fontsize=f_size)
        ax.spines['bottom'].set_gid("x_axis_line_" + column)
        ax.spines['left'].set_gid("y_axis_line_" + column)
        plt.rcParams['svg.fonttype'] = 'none'
        fig_name=f"replicate_1_vs_2_{s1}"
        plt.savefig(os.path.join(output_dir0, f"01_{column}_{fig_name}.svg"))
        plt.savefig(os.path.join(output_dir0, f"01_{column}_{fig_name}.png"))
        plt.close(fig)

correlation_df = pd.DataFrame({'Column_Name': column_names, 'Pearson_r': pearson_r_values})
# Calculate mean Pearson-R
mean_pearson_r = correlation_df['Pearson_r'].mean()
mean_row = pd.DataFrame({'Column_Name': ['Mean'], 'Pearson_r': [mean_pearson_r]})
correlation_df = pd.concat([correlation_df, mean_row], ignore_index=True)

# Export table to Excel
with pd.ExcelWriter(output_path) as writer:
    comparison_table.to_excel(writer, sheet_name=f"Replicate_comparison", index=False)
    correlation_df.to_excel(writer, sheet_name="Pearson R")
print(f"Tables saved in {output_path}.")

svg_paths = [os.path.join(output_dir0, f) for f in os.listdir(output_dir0) if f.endswith(f"replicate_1_vs_2_{s1}.svg")]
num_svgs = len(svg_paths)

def find_parent(root, element):
    """Finds the parent element of a given element in the XML tree."""
    for parent in root.iter():
        for child in parent:
            if child is element:
                return parent
    return None

def remove_axes_elements(root):
    """Removes axes-related g-elements based on the position in the grid."""
    elements_to_remove = []
    for elem in root.findall('.//{*}g'):
        group_id = elem.get('id')
        parent = find_parent(root, elem)
        parent_id=parent.get('id')
        if group_id:
            if "x_axis" in group_id and "matplotlib.axis_1" in parent_id:
                # Delete individual x-axis label
                elements_to_remove.append(elem)
            if "y_axis" in group_id and "matplotlib.axis_2" in parent_id:
                # Delete individual y-axis label
                elements_to_remove.append(elem)


    # Remove the collected elements after the iteration, using the parent node
    for elem_to_remove in elements_to_remove:
        parent = find_parent(root, elem_to_remove)
        if parent is not None:
            try:
                parent.remove(elem_to_remove)
            except ValueError:
                print(f"Warning: Element {elem_to_remove.get('id')} could not be found directly under the parent node.")
        else:
            print(f"Warning: Parent node of {elem_to_remove.get('id')} could not be found.")
    return root

if num_svgs == 0:
    print(f"No SVG files found in the folder: {output_dir0}")
else:
    svg_roots = []
    max_width = 0
    max_height = 0

    # 1. Parse all SVGs and get their dimensions, identify ax titles
    axes_elements = []
    parsed_svgs = []
    for path in svg_paths:
        try:
            tree = ET.parse(path)
            root = tree.getroot()
            parsed_svgs.append((root, path))

            width_str = root.get('width')
            height_str = root.get('height')
            width = 0
            height = 0

            if width_str and width_str.endswith('pt'):
                width = float(width_str[:-2])
            elif width_str:
                try:
                    width = float(width_str)
                except ValueError:
                    print(f"Warning: Invalid width format in {path}: {width_str}")

            if height_str and height_str.endswith('pt'):
                height = float(height_str[:-2])
            elif height_str:
                try:
                    height = float(height_str)
                except ValueError:
                    print(f"Warning: Invalid height format in {path}: {height_str}")

            max_width = max(max_width, width)
            max_height = max(max_height, height)

        except ET.ParseError as e:
            print(f"Error parsing {path}: {e}")

    if parsed_svgs:

        # total width and height with space for ax titles
        offset2= axes_offset*0.5
        total_width = cols * max_width + axes_offset - offset2 *(cols-1)
        total_height = rows * max_height + axes_offset #-2*offset2

        merged_svg = ET.Element('svg', {'width': str(total_width), 'height': str(total_height)})

        # common y-label
        ET.SubElement(merged_svg, 'text', {
            'x': str(axes_offset),
            'y': str(total_height / 2),
            'fill': 'black',
            'font-size': title_font_size,
            'text-anchor': 'middle',
            'transform': f'rotate(-90, {axes_offset *0.75}, {total_height / 2})'
        }).text = common_y_title

        # common x-label
        ET.SubElement(merged_svg, 'text', {
            'x': str(total_width / 2),
            'y': str(total_height - axes_offset / 4),
            'fill': 'black',
            'font-size': title_font_size,
            'text-anchor': 'middle'
        }).text = common_x_title

        # Position SVGs and remove axis elements
        for i, (root, path) in enumerate(parsed_svgs):
            col = i % cols
            row = i // cols

            # Modify the SVG file to remove axis elements
            modified_root = remove_axes_elements(root)

            if col > 8:
                x_offset_inner = col * max_width -offset2*(col-2)# Offset for the individual diagrams
            else:
                x_offset_inner = col * max_width +(axes_offset - (offset2 *col)) # Offset for the individual diagrams

            y_offset_inner = row * max_height

            group = ET.SubElement(merged_svg, 'g', {'transform': f'translate({x_offset_inner}, {y_offset_inner})'})

            for elem in modified_root:
                group.append(elem)

        merged_file_path = os.path.join(output_dir, output_filename)
        ET.ElementTree(merged_svg).write(merged_file_path)
        print(f"The merged SVG file with edited axes was saved under: {merged_file_path}")
    else:
        print("No valid SVG files found to merge.")