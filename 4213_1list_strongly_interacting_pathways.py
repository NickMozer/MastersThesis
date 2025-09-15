# Make a list of interacting pathways with a strong interaction above a 0.3 response fraction

import pandas as pd
import os

# Define the output directory based on the original script's logic
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)

# Load the comprehensive pathway interaction data
heatmap_data_all = pd.read_excel("04_pathway_interactions_cutoff_1.0.xlsx", "fraction_responses_over_1.0", index_col="pathway")

# Set the desired response fraction
response_fraction_cutoff = 0.3

# Initialize a list to store the pathway interactions that meet the criteria
interactions_at_cutoff = []

# Iterate through each pathway (row) in the heatmap data
for response_pathway, row in heatmap_data_all.iterrows():
    # Iterate through each knockdown pathway (column) and its response fraction
    for kd_pathway, response_fraction in row.items():
        # Check if the response fraction matches the desired cutoff
        if response_fraction >= response_fraction_cutoff:
            interactions_at_cutoff.append({
                "Response Pathway": response_pathway,
                "Knockdown Pathway": kd_pathway,
                "Response Fraction": response_fraction
            })

# Convert the list of interactions into a pandas DataFrame
df_interactions_at_cutoff = pd.DataFrame(interactions_at_cutoff)

# Define the output Excel file path
output_excel_filename = "13_pathway_interactions_at_0.3.xlsx"
output_excel_path = os.path.join(output_dir, output_excel_filename)

# Save the DataFrame to an Excel file
df_interactions_at_cutoff.to_excel(output_excel_path, sheet_name="Interactions_at_0.3", index=False)

print(f"Excel file '{output_excel_filename}' with pathway interactions at a response fraction of {response_fraction_cutoff} has been created in '{output_dir}'.")
