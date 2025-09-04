# Sort proteome data, calculate control means and errors, calculate fold changes of proteins upon knockdown with CRISPR from control means 

# Import libraries
import pandas as pd
import os

# Initialize output path
output_excel_path = "001_extract_proteome_control_and_crispri.xlsx"
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)
output_path = os.path.join(output_dir, output_excel_path)

# load dictionary tables from txt files
batch1_2_mapping = pd.read_csv("./batch1-234_may13/Ex-Ale-All.txt", sep="\t", usecols=["Name", "Experiment"])
batch4_mapping = pd.read_csv("./batch4-last50_nov14/EX-Ale-2-1Control.txt", sep="\t", usecols=["Name", "Experiment"])
samples_number_mapping=pd.read_excel("Samples in each batch.xlsx", "all_batches", usecols=["PROTEOME_ID", "STRAIN"])
samples_number_mapping["PROTEOME_ID"]=samples_number_mapping["PROTEOME_ID"].astype(str)

# Make Mapping-Dictionaries for Knockdown ID from txt
mapping_dict = pd.concat([batch1_2_mapping, batch4_mapping]).set_index("Name")["Experiment"].to_dict()
mapping_dict2 = samples_number_mapping.set_index("PROTEOME_ID")["STRAIN"].to_dict()

# Load proteomics data
excel_files = [
    ('./batch1-234_may13', "Ale-234.xlsx", "Ale-All_v2_PROTEIN"),
    ('./batch2-253_may17', "Ale-LastBatch.xlsx", "Ale-LastBatch_PROTEIN"),
    ('./batch4-last50_nov14', "Ale-2-1Control.xlsx", "Ale-2-1Control_PROTEIN")
]

# Sort data into CRISPR strains and control strains
tables = []
control_tables = []

for folder_path, excel_dataname, sheet_name in excel_files:
    file_path = os.path.join(folder_path, excel_dataname)
    table = pd.read_excel(file_path, sheet_name=sheet_name)
    # filter for columns containing "quantity" and "imp"
    selected_columns = ["proteinName", "geneName"] + [col for col in table.columns if "quantity" in col.lower()]
    filtered_table = table[selected_columns]
    selected_columns = ["proteinName", "geneName"] + [col for col in filtered_table.columns if "imp" in col.lower()]
    filtered_table = filtered_table[selected_columns]
    # transform data to the power of two
    filtered_table.iloc[:, 2:] = 2 ** filtered_table.iloc[:, 2:].astype(float)
    filtered_table.columns=filtered_table.columns.str.replace("Imp_Log2.", "")
    # Sort into batches, separate controls from CRISPRi strains
    if "batch4" in folder_path:
        # rename columns with dictionary
        filtered_table.rename(columns=mapping_dict, inplace=True)
        control_table = filtered_table.loc[:, ["proteinName", "geneName"] + [col for col in filtered_table.columns if "cntrl" in col.lower()]]
        filtered_table = filtered_table.loc[:, ["proteinName", "geneName"] + [col for col in filtered_table.columns if "cntrl" not in col.lower() and "name" not in col.lower()]]
        filtered_table.columns = filtered_table.columns.str.replace(r"^\d+_|_\d+$", "", regex=True)
    elif "batch1" in folder_path:
        # rename columns with dictionary
        filtered_table.rename(columns=mapping_dict, inplace=True)
        control_table = filtered_table.loc[:, ["proteinName", "geneName"] + [col for col in filtered_table.columns if "c" in col.lower()]]
        filtered_table = filtered_table.loc[:, ["proteinName", "geneName"] + [col for col in filtered_table.columns if "c" not in col.lower() and "name" not in col.lower()]]
        filtered_table.rename(columns=mapping_dict2, inplace=True)
    else:
        # rename columns with dictionary
        filtered_table.columns=filtered_table.columns.str.replace("Imp", "").str.replace("_", "").str.replace("Quantity", "").str.replace(".", "").str.replace("Log2", "").str.replace("CLib-", "").str.replace("raw", "")
        control_table = filtered_table.loc[:, ["proteinName", "geneName"] + [col for col in filtered_table.columns if "c" in col.lower()]]
        filtered_table = filtered_table.loc[:, ["proteinName", "geneName"] + [col for col in filtered_table.columns if "c" not in col.lower() and "name" not in col.lower()]]
        filtered_table.columns=filtered_table.columns.str.replace("1-", "").str.replace("re", "")
        filtered_table.rename(columns=mapping_dict2, inplace=True)
    control_tables.append(control_table)
    tables.append(filtered_table)

# Calculation of mean, SD, relSD of control strains per batch
for i, control_table in enumerate(control_tables):
    control_table[f'mean{i+1}'] = control_table.iloc[:, 2:].mean(axis=1)
    control_table[f'stdw{i+1}'] = control_table.iloc[:, 2:].std(axis=1)
    control_table[f'rsd{i+1}'] = (control_table[f'stdw{i+1}'] / control_table[f'mean{i+1}']) * 100

for i, filtered_table in enumerate(tables):
    filtered_table[f'mean{i+1}'] = control_tables[i][f'mean{i+1}']
    filtered_table[f'stdw{i+1}'] = control_tables[i][f'stdw{i+1}']
    filtered_table[f'rsd{i+1}'] = control_tables[i][f'rsd{i+1}']

# calculate fold changes
fc_tables = []
for i, fc_table in enumerate(tables):
    fc_table_fc = fc_table.copy()
    mean_col = f'mean{i+1}'  
    fc_table_fc.iloc[:, 2:-3] = fc_table_fc.iloc[:, 2:-3].div(fc_table_fc[mean_col], axis=0) 
    fc_tables.append(fc_table_fc)

# Export tables to Excel
with pd.ExcelWriter(output_path) as writer:
    for i, control_table in enumerate(control_tables):
        control_table.to_excel(writer, sheet_name=f"Control{i+1}", index=False)
    for i, filtered_table in enumerate(tables):
        filtered_table.to_excel(writer, sheet_name=f"Batch{i+1}", index=False)
    for i, fc_table in enumerate(fc_tables):
        fc_table.to_excel(writer, sheet_name=f"fc_Batch{i+1}", index=False)
print(f"Tables saved in {output_path}.")