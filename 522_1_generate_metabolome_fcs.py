# Create fold changes for metabolome data
# Merge pos and neg mode: Use all metabolites from negative mode, and then add those from pos mode. Note that everything with 0 was not detected by FI-MS
# Normalize to OD (intensity/OD)
# Create fold changes by dividing by the mean of controls
# Implement correct column names

# Import libraries
import pandas as pd
import numpy as np
import os

# Initialize output dir
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)

# Load dataframes
df_metab_neg = pd.read_excel("./mbx_neg_headers.xlsx", sheet_name="Sheet1")
df_metab_pos = pd.read_excel("./mbx_pos_headers.xlsx", sheet_name="Sheet1")
df_metab_ods = pd.read_excel("./Metabolome_samples.xlsx", sheet_name="Tabelle1")
df_ale_ods = pd.read_excel("./AA4_AleLib_growth_flask.xlsx", sheet_name="growth_flask_samples")
fc_table_prot = pd.read_excel("./005_log2_fcs_final.xlsx", sheet_name="Successfull_knockdowns")

# Select relevant columns (abbr and .mzxml columns)
selected_columns = ["abbr"] + [col for col in df_metab_neg.columns if ".mzxml" in col.lower()]
df_metab_neg_filter01 = df_metab_neg[selected_columns].copy()
df_metab_pos_filter = df_metab_pos[selected_columns].copy()
df_metab_neg_filter01.to_excel(os.path.join(output_dir, "210_metabolome_fc_table_neg_intens.xlsx"), index=False)
df_metab_neg_filter = df_metab_neg_filter01.copy()
print("Initial df_metab_neg_filter shape:", df_metab_neg_filter.shape)
print("Initial df_metab_pos_filter shape:", df_metab_pos_filter.shape)

# Merge pos and neg mode (replace zeros in neg with pos values)
# Iterate over columns that contain intensity data (ending with .mzxml)
mzxml_columns = [col for col in df_metab_neg_filter.columns if ".mzxml" in col.lower()]
mean_column = df_metab_neg_filter[mzxml_columns].mean(axis=1)
print("\n Step 1: Replacing zeros in neg mode with pos mode values")
# Identify rows in df_metab_neg_filter that are entirely zeros in mzxml_columns
# This creates a boolean Series: True if all mzxml columns in that row are zero.
all_zeros_mask = (df_metab_neg_filter[mzxml_columns] == 0).all(axis=1)

# Update 'abbr' column for these specific rows
# For the rows where all_zeros_mask is True, replace the 'abbr' from neg with 'abbr' from pos.
df_metab_neg_filter.loc[all_zeros_mask, 'abbr'] = df_metab_pos_filter.loc[all_zeros_mask, 'abbr']

# Use boolean mask for efficient replacement
# Find locations in df_metab_neg_filter where the value is 0
zero_mask = df_metab_neg_filter[mzxml_columns] == 0

# Replace these zeros with corresponding values from df_metab_pos_filter
df_metab_neg_filter[mzxml_columns] = df_metab_neg_filter[mzxml_columns].mask(zero_mask, df_metab_pos_filter[mzxml_columns])
df_metab_neg_filter.to_excel(os.path.join(output_dir, "211_metabolome_fc_table_filled_pos.xlsx"), index=False)
print("Completed replacement of zeros from negative mode with positive mode values.")
print("df_metab_neg_filter shape after replacement:", df_metab_neg_filter.shape)

# Remove all remaining rows with zeros in them
print("\n Step 2: Removing rows with any remaining zeros")
initial_rows = df_metab_neg_filter.shape[0]

# Check for zeros in the mzxml_columns and drop rows where any zero exists
# Use .any(axis=1) to check if any value in a row is 0 for the mzxml_columns
rows_to_keep = ~(df_metab_neg_filter[mzxml_columns] == 0).any(axis=1)
df_metab_filtered = df_metab_neg_filter[rows_to_keep].copy()
df_metab_filtered.to_excel(os.path.join(output_dir, "212_metabolome_fc_table_nozeros.xlsx"), index=False)
rows_removed = initial_rows - df_metab_filtered.shape[0]
print(f"Removed {rows_removed} rows containing zeros.")
print("df_metab_filtered shape after removing zeros:", df_metab_filtered.shape)

# Normalize intensity values by OD
print("\n Step 3: Normalizing intensity values by OD")

# Create a dictionary for quick OD lookup
od_mapping = df_metab_ods.set_index("Sample ID")["Sampling OD"].to_dict()

# Identify columns that need normalization (all mzxml columns)
cols_to_normalize = [col for col in df_metab_filtered.columns if ".mzxml" in col.lower()]

# Create a copy to store normalized data
df_metab_normalized = df_metab_filtered.copy()

# Normalize each relevant column
for col in cols_to_normalize:
    # The column name directly corresponds to the 'Sample' in df_metab_ods
    if col in od_mapping:
        od_value = od_mapping[col]
        if od_value != 0: # Avoid division by zero
            df_metab_normalized[col] = df_metab_normalized[col] / od_value
        else:
            print(f"Warning: OD value for {col} is zero. Skipping normalization for this column.")
            df_metab_normalized[col] = np.nan # Set to NaN if OD is zero to indicate problem
    else:
        print(f"Warning: OD value not found for {col}. Skipping normalization for this column.")
        df_metab_normalized[col] = np.nan # Set to NaN if OD is not found
df_metab_normalized.to_excel(os.path.join(output_dir, "213_metabolome_fc_table_normalized.xlsx"), index=False)
print("Completed normalization by OD values.")
print("df_metab_normalized shape:", df_metab_normalized.shape)

# Check for duplicate column names in df_fold_changes
print("\n Checking for duplicate column names in the final DataFrame")
strain_samples = {}
for index, row in df_ale_ods.iterrows():
    strain_name = str(row["strain"]) # Ensure strain_name is a string
    sample_name = row["worklist_met_code"]
    strain_samples[strain_name]=sample_name

print(len(strain_samples))
strains_to_keep = []
for strain in df_metab_normalized.columns:
    if strain == 'abbr':
        continue # Skip the abbr column
    parts = strain.split('_')
    strain_ID=parts[0]
    # Keep latest strains
    if strain_ID in strain_samples.values():
        strains_to_keep.append(strain)
    if strain.startswith("C"):
        strains_to_keep.append(strain)
print(len(strains_to_keep))
final_replicate_cols = ["abbr"]+[col for col in df_metab_normalized.columns if col in strains_to_keep]
df_single_replicates=df_metab_normalized[final_replicate_cols].copy()
df_single_replicates.to_excel(os.path.join(output_dir, "214_metabolome_fc_table_single_rep.xlsx"), index=False)

# Calculate fold changes by dividing by the mean of controls
print("\n Step 4: Calculating fold changes")

# Identify control columns and knockdown columns
# Control columns start with 'C' and end with '.mzxml'
cols_to_fc=[col for col in final_replicate_cols if "abbr" not in col]
control_cols = [col for col in cols_to_fc if col.startswith('C') and ".mzxml" in col.lower()]
# Knockdown columns are all other mzxml columns that are not controls
knockdown_cols = [col for col in cols_to_fc if col not in control_cols]

print(f"Found {len(control_cols)} control columns.")
print(f"Found {len(knockdown_cols)} knockdown columns.")

# Calculate the mean of control columns for each metabolite
# Ensure that all control columns are present in df_metab_normalized
# Also handle cases where a control column might have NaNs introduced by previous steps
control_means = df_single_replicates[control_cols].mean(axis=1)
control_stds = df_single_replicates[control_cols].std(axis=1)
control_rsds = (control_stds/ control_means) * 100

# Create a new DataFrame for fold changes
df_fold_changes = pd.DataFrame()
df_fold_changes["abbr"] = df_single_replicates["abbr"]

# Calculate fold changes for each knockdown column
for kd_col in knockdown_cols:
    # Avoid division by zero for control means
    # If control_means is 0 for a given metabolite, fold change will be NaN
    df_fold_changes[f"{kd_col.replace('.mzXML', '')}"] = df_metab_normalized[kd_col] / control_means.replace(0, np.nan)
df_fold_changes.columns=df_fold_changes.columns.str.replace("1", "").str.replace("2", "").str.replace("3", "").str.replace("4", "").str.replace("5", "").str.replace("6", "").str.replace("7", "").str.replace("8", "").str.replace("9", "").str.replace("0", "").str.replace("m_P", "").str.replace("-A", "")
df_fold_changes.columns=df_fold_changes.columns.str.replace("_P", "").str.replace("-P", "").str.replace("-B", "").str.replace("-C", "").str.replace("-D", "").str.replace("-E", "").str.replace("-F", "").str.replace("-G", "").str.replace("-H", "").str.replace("_pos", "").str.replace("atph", "atpH")
print("Completed fold change calculation.")
print("df_fold_changes shape:", df_fold_changes.shape)
print("First 5 rows of df_fold_changes:")
print(df_fold_changes.head())
df_fold_changes.to_excel(os.path.join(output_dir, "215_metabolome_fc_table.xlsx"), index=False)

# Delete dysfunctional columns
new_knockdown_cols= [col for col in df_fold_changes.columns if col not in control_cols]
kd_columns=["abbr"]+[col for col in fc_table_prot.columns if "gene" not in col and col in new_knockdown_cols]
df_func_kds=df_fold_changes.loc[:, kd_columns].copy()
df_func_kds.to_excel(os.path.join(output_dir, "216_metabolome_fc_table_functioning_cols.xlsx"), index=False)

# Calculate log2 fold changes for each knockdown column
df_fold_changes_log2=df_func_kds.copy()
df_fold_changes_log2.iloc[:,1:]=np.log2(df_fold_changes_log2.iloc[:,1:])

# Save the final fold change table to an Excel file
print(f"\nSaving fold changes to {output_dir}")
df_fold_changes_log2.to_excel(os.path.join(output_dir, "217_metabolome_fc_table_log2.xlsx"), index=False)

# Remove abbr suffixes 
suffixes_to_remove = r'\[M-H\]-|\[M\+H\]\+'
df_fc_log2_no_suffixes=df_fold_changes_log2.copy()
df_rsd=df_fc_log2_no_suffixes.copy()
df_fc_log2_no_suffixes['abbr'] = df_fc_log2_no_suffixes['abbr'].str.replace(suffixes_to_remove, '', regex=True)

# Keep low rsd metabolites only
df_rsd["mean_rsd"]=control_rsds
df_rsd['abbr'] = df_rsd['abbr'].str.replace(suffixes_to_remove, '', regex=True)
df_rsd.to_excel(os.path.join(output_dir, "218_metabolome_fc_table_rsds.xlsx"), index=False)
df_low_rsd=df_rsd.copy()
df_low_rsd = df_low_rsd[df_low_rsd["mean_rsd"] < 50]
df_low_rsd = df_low_rsd[[col for col in df_low_rsd.columns if "rsd" not in col]]
df_low_rsd.to_excel(os.path.join(output_dir, "219_metabolome_fc_table_low_rsds.xlsx"), index=False)

print("Analysis complete and results saved. Table 219 is the final result.")