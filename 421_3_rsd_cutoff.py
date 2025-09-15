# Filter for rsds below 50% to get proteins that can be measured well, identify outliers

# Import libraries
import pandas as pd
import os
import numpy as np

# cutoff at 50%
cutoff = 50

# Initialize output path
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_excel_path ="003_rsd_analysis.xlsx"
output_path = os.path.join(output_dir, output_excel_path)
output_excel_path2 = f"003_fcs_with_rsd_below_{cutoff}.xlsx"
output_path2 = os.path.join(output_dir, output_excel_path2)

# Load rsd table
rsd_table = pd.read_excel("./002_log2_fcs_and_errors.xlsx", "just_rsds")

# Filter rows for mean_rsd or single rsd, cutoff at 50%
outlier_table = rsd_table[rsd_table["mean_rsd"] > cutoff]
rsd_mean_below50_table = rsd_table[rsd_table["mean_rsd"] < cutoff]
rsd_all_below50_table= rsd_mean_below50_table[rsd_mean_below50_table["rsd1"] < cutoff]
rsd_all_below50_table= rsd_all_below50_table[rsd_all_below50_table["rsd2"] < cutoff]
rsd_all_below50_table= rsd_all_below50_table[rsd_all_below50_table["rsd3"] < cutoff]

# Make outlier table
rsd_missing_rows = rsd_mean_below50_table.merge(
    rsd_all_below50_table, 
    on=list(rsd_mean_below50_table.columns),
    how="left",
    indicator=True
).query('_merge == "left_only"').drop(columns=["_merge"])

# Save tables to excel
with pd.ExcelWriter(output_path) as writer:      
    outlier_table.to_excel(writer, sheet_name=f"rsd_outliers", index=False)
    rsd_mean_below50_table.to_excel(writer, sheet_name=f"rsd_mean_below{cutoff}",index=False)
    rsd_all_below50_table.to_excel(writer, sheet_name=f"rsd_all_below{cutoff}",index=False)
    rsd_missing_rows.to_excel(writer, sheet_name=f"rsd_other_outliers",index=False)
print(f" Rsd tables saved in {output_path}.")

# Load fc, rsd tables
fc_table = pd.read_excel("./002_fcs_and_errors.xlsx", "just_fcs")

# Merge rsd_table and fc_table on geneName, remove rsds
low_error_table = pd.merge(fc_table, rsd_all_below50_table, on='geneName')
low_error_table = low_error_table.loc[:, [col for col in low_error_table.columns if "rsd" not in col.lower()]]

# Reshape the data for boxplot
fc_values = low_error_table.iloc[:, 1:].values.flatten()

# Log2 only fc-table
fcs_table_log2=low_error_table.copy()
fcs_table_log2.iloc[:,1:]=np.log2(fcs_table_log2.iloc[:,1:])

# Save table to Excel
with pd.ExcelWriter(output_path2) as writer:      
    low_error_table.to_excel(writer, sheet_name=f"fcs_with_rsd_below_{cutoff}", index=False)
    fcs_table_log2.to_excel(writer, sheet_name=f"log2_fcs_with_rsd_below_{cutoff}", index=False)
print(f"Fc tables saved in {output_path2}.")
