# Create a dendrogram to cluster the correlation of all dataset 2 proteomes

import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram, set_link_color_palette
from scipy.spatial.distance import pdist

input_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(input_dir)

output_dir_clustering = os.path.join(input_dir, "proteome_clusters_prot")
os.makedirs(output_dir_clustering, exist_ok=True) 

fc_table_prot = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")

kd_columns_prot = [col for col in fc_table_prot.columns if "geneName" not in col]

# Prepare Data for Clustering
# Set 'geneName' as the DataFrame index. 
df_clustering_prot = fc_table_prot.set_index('geneName')[kd_columns_prot]

# Handle missing values (NaNs).
initial_nan_count_prot = df_clustering_prot.isnull().sum().sum()
if initial_nan_count_prot > 0:
    print(f"Warning: {initial_nan_count_prot} NaN values found in the data. Filling them with 0.")
    df_clustering_prot = df_clustering_prot.fillna(0)
else:
    print("No NaN values found in the data. Data is clean.")

# Transpose the DataFrame for clustering.
df_clustering_prot_transposed = df_clustering_prot.T
print(f"Data prepared for clustering. Shape: {df_clustering_prot_transposed.shape} (Knockdowns x Proteins)")

# Perform Hierarchical Clustering and Visualize Dendrogram
print("Performing hierarchical clustering and generating dendrogram...")

# Calculate the pairwise distances_prot between knockdowns.
# 'correlation' metric calculates 1 - Pearson correlation coefficient. This is the distance metric, where 0 means perfect correlation and 2 means perfect anti-correlation.
distances_prot = pdist(df_clustering_prot_transposed, metric='correlation')

# Perform hierarchical clustering using the 'average' linkage method.
# This method calculates the distance between two clusters_prot as the average distance between all observations in the two clusters_prot.
row_linkage_prot = linkage(distances_prot, method='average')

# Initialize a DataFrame to store the knockdown names and their assigned clusters_prot for each threshold.
cluster_assignments_prot = pd.DataFrame({'Knockdown': df_clustering_prot_transposed.index})

# Define the desired Pearson correlation thresholds.
thresholds = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
custom_colors = ("darkgreen", "darkviolet", "firebrick", "mediumblue", "chocolate")

set_link_color_palette(list(custom_colors))
# Loop through each threshold to determine and store cluster assignments.
plt.rcParams['font.size'] = 24
for threshold in thresholds:
    fig_width = 21 # Desired figure width remains 21
    fig_height = 30 # Desired figure height remains 30

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Adjust subplot parameters to make space for labels on the right
    fig.subplots_adjust(left=0.00, right=0.85, top=1.0, bottom=0.00) 
    plt.xlabel("Distance (1 - Pearson Correlation)")
    plt.ylabel("CRISPRi Knockdowns")

    min_pearson_correlation = threshold
    distance_threshold = 1 - min_pearson_correlation

    current_clusters = fcluster(row_linkage_prot, distance_threshold, criterion='distance')
    
    # Create a mapping from original data index to its cluster ID
    original_idx_to_cluster_id = {i: cluster_id for i, cluster_id in enumerate(current_clusters)}

    # Plot the dendrogram WITHOUT labels
    dendro = dendrogram(
        row_linkage_prot,
        labels=df_clustering_prot_transposed.index.tolist(), # Still pass labels for internal use
        orientation='left',
        no_labels=True, # Suppress dendrogram's own labels
        color_threshold=distance_threshold,
        above_threshold_color='black',
        ax=ax
    )

    # Manually add and adjust text labels
    texts = []
    # Get the indices of the leaves in the order they appear on the dendrogram
    leaf_indices = dendro['leaves']
    # Get the original labels based on these indices
    original_labels = df_clustering_prot_transposed.index.tolist()
    sorted_labels = [original_labels[i] for i in leaf_indices]

    # The y-coordinates of the leaf labels are typically evenly spaced by dendrogram
    # The x-coordinate for labels will be just to the right of the dendrogram's extent.
    max_x_coord = ax.get_xlim()[1]  # Get the current max x-limit of the plot

    # Calculate y-coordinates for labels
    num_leaves = len(df_clustering_prot_transposed)
    y_coords = [10 * i + 5 for i in range(num_leaves)] # Default spacing for leaves in scipy dendrogram

    # Map the sorted_labels to their actual y-coordinates on the plot.
    # Loop through the leaves in the order they are plotted by dendrogram
    counter=1
    for i, leaf_idx in enumerate(dendro['leaves']):
        label = original_labels[leaf_idx]
        
        # Get the cluster ID for this specific leaf (original data point)
        cluster_id = original_idx_to_cluster_id[leaf_idx]
        # x-coordinate for the label: just to the right of the dendrogram's farthest point
        # A small offset is added for visual separation
        if counter==1:
            counter=2
            x_pos = max_x_coord - 0.12
        elif counter==2:
            counter=3
            x_pos = max_x_coord - 0.06
        elif counter==3:
            counter=1
            x_pos = max_x_coord
        else:
            x_pos = max_x_coord
        
        temp_fig, temp_ax = plt.subplots(figsize=(1,1)) # Small, temporary figure
        temp_dendro = dendrogram(
            row_linkage_prot,
            orientation='left',
            color_threshold=distance_threshold,
            above_threshold_color='black',
            no_plot=True, # Do not plot to screen
            labels=df_clustering_prot_transposed.index.tolist() # Need labels to get leaves_color_list
        )
        plt.close(temp_fig) # Close the temporary figure immediately

        original_idx_to_leaf_color = {
            temp_dendro['leaves'][k]: temp_dendro['leaves_color_list'][k]
            for k in range(len(temp_dendro['leaves']))
        }
        # Get the color for the current leaf in the main dendrogram plot
        label_color = original_idx_to_leaf_color[leaf_idx]
        # y-coordinate for the label: based on its position in the dendrogram's leaf order
        # The y-positions are generated by scipy; for 'left' orientation, they are typically
        y_pos = 10 * i + 5 

        texts.append(plt.text(x_pos, y_pos, label, va='center', ha='left', color=label_color))
    for pos in ['right', 'top']:
        plt.gca().spines[pos].set_visible(False)

    plt.tight_layout()

    # Save the generated dendrogram to a file.
    dendrogram_output_path1 = os.path.join(output_dir_clustering, f"knockdown_proteome_dendrogram_{min_pearson_correlation}.svg")
    dendrogram_output_path2 = os.path.join(output_dir_clustering, f"knockdown_proteome_dendrogram_{min_pearson_correlation}.png")
    plt.savefig(dendrogram_output_path1, dpi=300, bbox_inches='tight')
    plt.savefig(dendrogram_output_path2, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Dendrogram saved to: {dendrogram_output_path1}")

    # Extract Cluster Assignments for Multiple Thresholds
    print("Extracting cluster assignments for multiple correlation thresholds...")

    clusters_prot = fcluster(row_linkage_prot, distance_threshold, criterion='distance')
    cluster_assignments_prot[f'Cluster_Corr_{min_pearson_correlation}'] = clusters_prot
    num_clusters_prot_formed = len(set(clusters_prot))
    print(f"   - For minimum Pearson correlation {min_pearson_correlation}: {num_clusters_prot_formed} clusters_prot formed.")

# Sort the assignments by the last cluster column
cluster_assignments_prot = cluster_assignments_prot.sort_values(by=cluster_assignments_prot.columns[9]).reset_index(drop=True)

# Save the cluster assignments to an Excel file.
cluster_assignments_prot_output_path = os.path.join(output_dir_clustering, "knockdown_cluster_assignments_prot_multi_threshold.xlsx")
cluster_assignments_prot.to_excel(cluster_assignments_prot_output_path, index=False)

print(f"Cluster assignments for multiple thresholds saved to: {cluster_assignments_prot_output_path}")
