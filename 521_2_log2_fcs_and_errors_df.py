# Sort fc data alphabetically, merge batches to one df, make one df without mean columns, make a log2 fc df, make dfs with just means or just rsds

# Import libraries
import pandas as pd
import numpy as np
import os

# Initialize output path
output_excel_path = "002_log2_fcs_and_errors.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_path = os.path.join(output_dir, output_excel_path)

# Load fc-data
excel_files = [
    (output_dir, f"001_extract_proteome_control_and_crispri.xlsx", "fc_Batch1"),
    (output_dir, f"001_extract_proteome_control_and_crispri.xlsx", "fc_Batch2"),
    (output_dir, f"001_extract_proteome_control_and_crispri.xlsx", "fc_Batch3")
]
# Sort CRISPR strains alphabetically
fc_tables = []
for folder_path, excel_dataname, sheet_name in excel_files:
    file_path = os.path.join(folder_path, excel_dataname)
    fc_table = pd.read_excel(file_path, sheet_name=sheet_name)
    fc_table = fc_table[
        list(fc_table.columns[:2]) + sorted(fc_table.columns[2:-3]) + list(fc_table.columns[-3:])
    ]
    fc_tables.append(fc_table)

# Merge fc Tables to a new one, double columns get _01,_02, _03; proteinName column eliminated
fc_one_table = pd.merge(fc_tables[0], fc_tables[1], on='geneName', how='outer', suffixes=('_01', '_02'))
fc_one_table = pd.merge(fc_one_table, fc_tables[2], on='geneName', how='outer', suffixes=('_01', '_03'))
fc_one_table=fc_one_table.loc[:, [col for col in fc_one_table.columns if "protein" not in col.lower()]]
fc_one_table=fc_one_table.sort_values(by='geneName')

# Only fcs table without mean, stdw, rsd
only_fcs_table=fc_one_table.loc[:, [col for col in fc_one_table.columns if "mean" not in col.lower()]]
only_fcs_table=only_fcs_table.loc[:, [col for col in only_fcs_table.columns if "rsd" not in col.lower()]]
only_fcs_table=only_fcs_table.loc[:, [col for col in only_fcs_table.columns if "stdw" not in col.lower()]]
only_fcs_table = only_fcs_table[
    list(only_fcs_table.columns[:1]) + sorted(only_fcs_table.columns[1:])
]
# Log2 only fc-table
only_fcs_table_log2=only_fcs_table.copy()
only_fcs_table_log2.iloc[:,1:]=np.log2(only_fcs_table_log2.iloc[:,1:])

# Only rsds table
only_rsds_table=fc_one_table.loc[:, ["geneName"] + [col for col in fc_one_table.columns if "rsd" in col.lower()]].copy()
only_rsds_table[f'mean_rsd'] = only_rsds_table[['rsd1', 'rsd2', 'rsd3']].mean(axis=1)

def median_ignore_nan(row):
    """Calculate Median of a row and ignore NaN-values."""
    valid_values = row.dropna().values
    if len(valid_values) < 3:
        return np.nan
    return np.median(valid_values)

only_rsds_table['median_rsd'] = only_rsds_table[['rsd1', 'rsd2', 'rsd3']].apply(median_ignore_nan, axis=1)

# Only means table
only_means_table=fc_one_table.loc[:, ["geneName"] + [col for col in fc_one_table.columns if "mean" in col.lower()]].copy()
only_means_table[f'mean_mean'] = only_means_table[['mean1', 'mean2', 'mean3']].mean(axis=1)
only_means_table[f'median_mean'] = only_means_table[['mean1', 'mean2', 'mean3']].apply(median_ignore_nan, axis=1)

# Export tables to Excel
with pd.ExcelWriter(output_path) as writer:      
    fc_one_table.to_excel(writer, sheet_name=f"all_fcs",index=False)
    only_fcs_table.to_excel(writer, sheet_name=f"just_fcs",index=False)
    only_fcs_table_log2.to_excel(writer, sheet_name=f"log2_fcs",index=False)
    only_means_table.to_excel(writer, sheet_name=f"just_means",index=False)
    only_rsds_table.to_excel(writer, sheet_name=f"just_rsds",index=False)
print(f"Tables saved in {output_path}.")