# Violin plots of all 281 knockdown proteomes with automatic distribution classification

# Import libraries
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats 
import numpy as np

# Initialize date and output path
input_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(input_dir)

# Define output directories
output_dir_plots = os.path.join(input_dir, "kd_violins_classified")
output_dir_classifications = input_dir
os.makedirs(output_dir_plots, exist_ok=True)

# Load fc_table 
fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")

# Identify knockdown columns (excluding 'geneName' identifier column)
kd_columns = [col for col in fc_table.columns if "geneName" not in col and "Unnamed" not in col]

# Dictionary to store classifications for each knockdown
distribution_classifications = {}

# Iterate through each knockdown proteome
for kd in kd_columns:
    # Extract fold change values for the current knockdown
    fc_values = fc_table.loc[:, kd].values.flatten()
    
    # Remove NaN values, as most statistical tests do not handle them
    fc_values = fc_values[~np.isnan(fc_values)]

    # Skip if there's insufficient data for statistical analysis
    if len(fc_values) < 3: # Shapiro-Wilk test requires at least 3 data points
        print(f"Skipping statistical analysis for '{kd}' due to insufficient data points ({len(fc_values)}).")
        distribution_classifications[kd] = "Insufficient Data"
        # Still plot if there's any data, but without classification in title
        if len(fc_values) == 0:
            print(f"No data to plot for '{kd}'. Skipping plot generation.")
            continue
    else:
        # Calculate statistical properties
        skewness = stats.skew(fc_values) # Measures the asymmetry of the distribution
        kurtosis = stats.kurtosis(fc_values) # Measures the 'tailedness' of the distribution (Fisher's definition, excess kurtosis)
        
        # Perform Shapiro-Wilk test for normality
        # Null hypothesis: The data is drawn from a normal distribution. A small p-value (e.g., < 0.05) indicates rejection of normality.
        shapiro_statistic, shapiro_p_value = stats.shapiro(fc_values)

        # Classify the distribution based on statistical properties
        dist_type = "Unclassified" # Default type

        # Define thresholds for classification
        SKEW_THRESHOLD_MODERATE = 0.5
        SKEW_THRESHOLD_HIGH = 1.0
        KURT_THRESHOLD = 1.0 # For excess kurtosis

        # Check for approximate normality first
        if shapiro_p_value > 0.05 and \
           abs(skewness) < SKEW_THRESHOLD_MODERATE and \
           abs(kurtosis) < KURT_THRESHOLD:
            dist_type = "Approximately Normal"
        else:
            # If not normal, describe skewness and kurtosis characteristics
            skew_description = ""
            if skewness >= SKEW_THRESHOLD_HIGH:
                skew_description = "Highly Right-skewed"
            elif skewness >= SKEW_THRESHOLD_MODERATE:
                skew_description = "Moderately Right-skewed"
            elif skewness <= -SKEW_THRESHOLD_HIGH:
                skew_description = "Highly Left-skewed"
            elif skewness <= -SKEW_THRESHOLD_MODERATE:
                skew_description = "Moderately Left-skewed"
            
            kurt_description = ""
            if kurtosis >= KURT_THRESHOLD:
                kurt_description = "Leptokurtic (Peaked)" # More peaked than normal
            elif kurtosis <= -KURT_THRESHOLD:
                kurt_description = "Platykurtic (Flat)" # Flatter than normal
            
            # Combine descriptions
            if skew_description and kurt_description:
                dist_type = f"{skew_description}, {kurt_description}"
            elif skew_description:
                dist_type = skew_description
            elif kurt_description:
                dist_type = kurt_description
            else:
                # If it failed normality but didn't fit simple skew/kurtosis thresholds
                dist_type = "Complex/Non-Normal" 
        
        distribution_classifications[kd] = dist_type

    # Plotting code (modified to include classification)
    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_subplot(1, 1, 1)

    # Violin width relative to figure width
    violin_width = 0.8
    ax.set_xticks([0])
    ax.set_xlim(-0.75, 0.75) # Space left and right

    # Violin Plot for non-pathway related responses (gray)
    sns.violinplot(data=fc_values,
                   color='darkgray', alpha=0.5, linewidth=0.5, zorder=1, width=violin_width, ax=ax)
    plt.xlabel(f'{kd}\n(n={len(fc_values)})')
    plt.ylabel('Log2 fold change')

    # Add the classification to the plot title
    if kd in distribution_classifications:
        plt.title(f"Distribution Type: {distribution_classifications[kd]}")
    else:
        plt.title(f"{kd} Distribution")

    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    plt.tick_params(axis='y', which='both', right=False, left=True, labelleft=True)
    for pos in ['right', 'top']:
        plt.gca().spines[pos].set_visible(False)
    plt.xticks(ha='right')

    plt.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'

    # Save the violin plot with classification in the filename
    if "High" and "Left" in dist_type:
        output_path_plot = os.path.join(output_dir_plots, f"./highleft/27_proteome_{kd}_violin_classified.png")
    elif "High" and "Right" in dist_type:
        output_path_plot = os.path.join(output_dir_plots, f"./highright/27_proteome_{kd}_violin_classified.png")
    elif "Mod" and "Left" in dist_type:
        output_path_plot = os.path.join(output_dir_plots, f"./modleft/27_proteome_{kd}_violin_classified.png")
    elif "Mod" and "Right" in dist_type:
        output_path_plot = os.path.join(output_dir_plots, f"./modright/27_proteome_{kd}_violin_classified.png")
    else:
        output_path_plot = os.path.join(output_dir_plots, f"./normal/27_proteome_{kd}_violin_classified.png")
    plt.savefig(output_path_plot)
    plt.close()

    print(f"Violin plot for '{kd}' (Type: {distribution_classifications.get(kd, 'N/A')}) saved to {output_path_plot}")

# Save all distribution classifications to a xlsx file
classification_df = pd.DataFrame.from_dict(distribution_classifications, orient='index', columns=['Distribution Type'])
classification_df.index.name = 'Knockdown'
output_path_xlsx = os.path.join(output_dir_classifications, "proteome_distribution_types.xlsx")
classification_df.to_excel(output_path_xlsx)

print(f"\nAll distribution classifications summarized in: {output_path_xlsx}")
print("\nProcess complete.")
