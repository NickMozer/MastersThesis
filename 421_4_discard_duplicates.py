# Get rid of duplicate CRISPRi columns by choosing the one with the lower self-response log2 fc or else one measured first

#Import libraries
import pandas as pd
import os

# Initialize output path
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
input_excel_path = f"003_fcs_with_rsd_below_50.xlsx"
output_excel_path = f"004_merged_crispri_log2_fcs.xlsx"
output_path = os.path.join(output_dir, output_excel_path)

# Load fc table
log2fc_table = pd.read_excel(input_excel_path, sheet_name="log2_fcs_with_rsd_below_50")

# Identify Crispri columns with number suffix
crispri_cols = [col for col in log2fc_table.columns if col != 'geneName' and '_' in col and col.split('_')[-1].isdigit()]

# Create a dictionary to store the base names of the CrispR strains
crispri_base_names = set([col.split('_')[0] for col in crispri_cols])

# Iterate through the base names of the Crispri strains
for base_name in crispri_base_names:
    # Find all columns that belong to this Crispri strain
    duplicate_cols = [col for col in crispri_cols if col.startswith(f"{base_name}_")]

    if len(duplicate_cols) > 1:
        print(f"Work on duplicates of: {base_name}")

        # Find the rows in which the ‘geneName’ column matches the base name
        gene_name_rows = log2fc_table[log2fc_table['geneName'] == base_name].index

        if not gene_name_rows.empty:
            # Only consider the first row found
            index_to_compare = gene_name_rows[0]

            fold_changes = {}
            for col in duplicate_cols:
                if col in log2fc_table.columns and not pd.isna(log2fc_table.at[index_to_compare, col]):
                    fold_changes[col] = log2fc_table.at[index_to_compare, col]

            if fold_changes:
                # Find the column with the lower fold change in this row
                best_col = min(fold_changes, key=fold_changes.get)
                # Transfer the value of the best column to a new column without suffix
                log2fc_table[base_name] = log2fc_table[best_col]
            else:
                # Fallback: If there are no fold change values in the duplicates, keep the column with the lower suffix
                best_col = min(duplicate_cols, key=lambda x: int(x.split('_')[-1]))
                log2fc_table[base_name] = log2fc_table[best_col]

            # Discard the original duplicate columns
            log2fc_table.drop(columns=duplicate_cols, inplace=True, errors='ignore')

            # Update the list of Crispri columns
            crispri_cols = [col for col in log2fc_table.columns if col != 'geneName' and '_' in col and col.split('_')[-1].isdigit()]

        else:
            print(f"No row with 'geneName' == '{base_name}' found for comparison of duplicates.")
            # Fallback: If there are no fold change values in the duplicates, keep the column with the lower suffix
            best_col = min(duplicate_cols, key=lambda x: int(x.split('_')[-1]))
            log2fc_table[base_name] = log2fc_table[best_col]

            # Discard the original duplicate columns
            log2fc_table.drop(columns=duplicate_cols, inplace=True, errors='ignore')
            # Update the list of Crispri columns
            crispri_cols = [col for col in log2fc_table.columns if col != 'geneName' and '_' in col and col.split('_')[-1].isdigit()]
    elif len(duplicate_cols) == 1:
        # If only one duplicate exists, rename it and remove the suffix
        log2fc_table.rename(columns={duplicate_cols[0]: base_name}, inplace=True, errors='ignore')
        crispri_cols = [col for col in log2fc_table.columns if col != 'geneName' and '_' in col and col.split('_')[-1].isdigit()]

log2fc_table = log2fc_table[list(log2fc_table.columns[:1]) + sorted(log2fc_table.columns[1:])]

# Save the cleaned table in a new Excel file
with pd.ExcelWriter(output_path) as writer:
    log2fc_table.to_excel(writer, sheet_name=f"log2_fcs_merged_crispri", index=False)

print(f"Fc tables saved in {output_path}.")
