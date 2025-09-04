# Filter out dysfunctional targets from original fc table
# Condition 1: no positive self-responses to knockdowns
# Condition 2: no OD below 0.3

# Import libraries
import pandas as pd
import os

# Define output directory
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_excel_path = f"005_log2_fcs_final.xlsx"
output_path = os.path.join(output_dir, output_excel_path)

# Load the fold change tables
fc_table = pd.read_excel("./004_merged_crispri_log2_fcs.xlsx", "log2_fcs_merged_crispri")
od_table = pd.read_excel("./AA4_AleLib_growth_flask.xlsx", "growth_flask_samples")

knockdown_genes = [col for col in fc_table.columns if col != 'geneName']
geneName_column=fc_table["geneName"]

knockdown_response = {}

# Apply Condition 1: no positive self-responses to knockdowns
for gene in knockdown_genes:
    response_found = False
    for index, name in geneName_column.items():
        if gene == name:
            knockdown_response[gene] = fc_table.loc[index, gene]
            response_found = True
            break
        elif gene == "trpD" and name == "trpGD":
            knockdown_response[gene] = fc_table.loc[index, gene]
            response_found = True
            break
    if not response_found:
        knockdown_response[gene] = -float('inf') # Treat as successful if no direct response measured

knockdown_response_df = pd.DataFrame(list(knockdown_response.items()), columns=['Knockdown', 'log2FoldChange'])

successfull_knockdowns_cols = ['geneName']
dysfunctional_knockdowns_cols = ['geneName']

for gene in knockdown_genes:
    if knockdown_response[gene] <= 0:
        successfull_knockdowns_cols.append(gene)
    else:
        dysfunctional_knockdowns_cols.append(gene)

Successfull_knockdowns = fc_table[successfull_knockdowns_cols]
Dysfunctional_knockdowns = fc_table[dysfunctional_knockdowns_cols]

# Apply Condition 2: No OD below 0.3

# Prepare a dictionary to store OD values for each target
target_ods = {}
numbers=(1,2,3,4,5,6,7,8,9,0)
for index, row in od_table.iterrows():
    sample_name = str(row["strain"]) # Ensure sample_name is a string
    od_value = row["final_OD"]

    # Check if it's a control strain (starts with 'C') or if OD is missing/invalid
    if not sample_name.startswith('C') and pd.notna(od_value):
        # Extract the target name by splitting at the first underscore
        parts = sample_name.split('_')
        if parts:
            target_name = parts[0]
            for number in numbers:
                target_name=target_name.replace(f"{number}m", "")
                target_name=target_name.replace(f"{number}", "")
            target_ods[target_name] = od_value # second Sample counts for Samples with two replicates

# Identify genes from 'Successfull_knockdowns' that have an OD below 0.3
genes_to_remove_due_to_od = []
for gene in Successfull_knockdowns.columns:
    if gene == 'geneName':
        continue # Skip the geneName column

    # Check if the gene is in our processed ODs and if its OD is below 0.3
    if gene in target_ods and target_ods[gene] < 0.3:
        genes_to_remove_due_to_od.append(gene)

# Filter 'Successfull_knockdowns' to remove columns with low OD
# These genes move from 'Successfull_knockdowns' to 'Dysfunctional_knockdowns'
final_successfull_knockdowns_cols = [col for col in successfull_knockdowns_cols if col not in genes_to_remove_due_to_od]
Final_Successfull_knockdowns = Successfull_knockdowns[final_successfull_knockdowns_cols]
Final_Successfull_knockdowns.columns=Final_Successfull_knockdowns.columns.str.replace("1", "")

# Add the genes removed due to low OD to 'Dysfunctional_knockdowns'
# Ensure no duplicates if a gene was already in dysfunctional_knockdowns_cols from condition 1
for gene in genes_to_remove_due_to_od:
    if gene not in dysfunctional_knockdowns_cols:
        dysfunctional_knockdowns_cols.append(gene)

# Recreate Dysfunctional_knockdowns with the newly identified dysfunctional targets
Dysfunctional_knockdowns_final = fc_table[dysfunctional_knockdowns_cols]

# Export tables to Excel
with pd.ExcelWriter(output_path) as writer:
    knockdown_response_df.to_excel(writer, sheet_name=f"Knockdown_own_response", index=False)
    Final_Successfull_knockdowns.to_excel(writer, sheet_name=f"Successfull_knockdowns", index=False)
    Dysfunctional_knockdowns_final.to_excel(writer, sheet_name=f"Dysfunctional_knockdowns", index=False)
print(f"Tables saved in {output_path}.")