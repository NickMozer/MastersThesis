# Correlations of all 281 knockdown proteomes 

# Import libraries
import pandas as pd
import os
from scipy.stats import pearsonr

# Initialize date and output path
input_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(input_dir)
output_dir_classifications = input_dir

# Load fc, operon tables
operon_table = pd.read_excel("./operons_ecoli.xlsx", "Tabelle1")
filtered_operon_table = operon_table[operon_table['CsiR'].astype(str).str.contains('-')].copy()

pathway_table = pd.read_excel("./20250327_pathways_reactants_ecoli.xlsx", "subpathways")
pathway_table = pathway_table.rename(columns={'GeneName': 'geneName', 'Pathway': 'pathway'})
ecocyc_gene_pathways = pathway_table.groupby('geneName')['pathway'].apply(list).to_dict()

regulon_table = pd.read_excel("./Regulons_regulondb.xlsx", "Regulator_gene")
regulon_table = regulon_table.rename(columns={'5)regulatedName': 'geneName', '3)RegulatorGeneName': 'regulon'})
ecocyc_gene_regulons = regulon_table.groupby('geneName')['regulon'].apply(list).to_dict()

fc_table = pd.read_excel("./005_log2_fcs_final.xlsx", "Successfull_knockdowns")
kd_columns=[col for col in fc_table.columns if "geneName" not in col]
fc_values_general=fc_table.iloc[:, 1:].values.flatten()

distribution_classifications = {}
relation_n=0
relation_o=0
relation_p=0
relation_r=0

modneg=0
noco =0
veweco=0
vewenegco=0
weco=0
wenegco=0
modco=0

# function to test wether two genes are in the same operon
def are_genes_in_same_operon(gene1, gene2):
    #if gene1 in filtered_operon_table['CsiR'].split('-') and gene2 in filtered_operon_table['CsiR'].split('-'):
    for _, row in filtered_operon_table.iterrows():
        if isinstance(row['CsiR'], str) and '-' in row['CsiR']:
            operon_genes = row['CsiR'].split('-')
            if gene1 in operon_genes and gene2 in operon_genes:
                return True
    return False

def are_genes_in_same_pathway(gene1, gene2):
    if gene1 in ecocyc_gene_pathways and gene2 in ecocyc_gene_pathways:
        ecocyc_pathways_gene1 = set(ecocyc_gene_pathways[gene1])
        ecocyc_pathways_gene2 = set(ecocyc_gene_pathways[gene2])
        if ecocyc_pathways_gene1.intersection(ecocyc_pathways_gene2):
            return True
    return False

def are_genes_in_same_regulon(gene1, gene2):
    if gene1 in ecocyc_gene_regulons and gene2 in ecocyc_gene_regulons:
        ecocyc_regulons_gene1 = set(ecocyc_gene_regulons[gene1])
        ecocyc_regulons_gene2 = set(ecocyc_gene_regulons[gene2])
        if ecocyc_regulons_gene1.intersection(ecocyc_regulons_gene2):
            return True
    return False

for i, kd1 in enumerate(kd_columns):
    kd_columns2_slice = kd_columns[i+1:]
    
    for kd2 in kd_columns2_slice:

        # Reshape the data for boxplot
        distribution_classifications[f"{kd1}_{kd2}"] = ("Insufficient Data", "NoType", "NoCommons")
        fc_values_kd1 = fc_table.loc[:, kd1].values.flatten()
        fc_values_kd2 = fc_table.loc[:, kd2].values.flatten()

        # Pearson R
        correlation, _ = pearsonr(fc_values_kd1, fc_values_kd2)

        if abs(correlation) < 0.01:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "No_corr", "")
                noco+=1
                continue
        elif abs(correlation) < 0.2:
            if correlation > 0:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Very_weak_corr", "")
                veweco+=1
                continue
            else:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Very_weak_neg_corr", "")
                vewenegco+=1
                continue
        elif abs(correlation) < 0.4:
            if correlation > 0:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Weak_corr", "")
                weco+=1
                continue
            else:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Weak_neg_corr", "")
                wenegco+=1
                continue
        elif abs(correlation) < 0.6:
            if correlation > 0:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Moderate_corr", "")
                modco+=1
                continue
            else:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Moderate_neg_corr", "")
                modneg+=1
                continue 
        if are_genes_in_same_operon(kd1, kd2):
            relation_o+=1
            relation="Same Operon"
        elif are_genes_in_same_pathway(kd1, kd2):
            relation_p+=1
            relation="Same pathway"
        elif are_genes_in_same_regulon(kd1, kd2):
            relation_r+=1
            relation="Same regulon"
        else:
            relation_n+=1
            relation="No Relation"

        if abs(correlation) < 0.8:
            if correlation > 0:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Strong_corr", relation)
            else:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Strong_neg_corr", relation)
        elif abs(correlation) <= 1:
            if correlation > 0:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Very_strong_corr", relation)
            else:
                distribution_classifications[f"{kd1}_{kd2}"] = (correlation, "Very_strong_neg_corr", relation)
        else:
            print("Error: Correlation above 1")

    print(f"Did all corrs for {kd1}")

print("Number of correlations:")
print(len(distribution_classifications))
# Save all distribution classifications to a xlsx file
classification_df = pd.DataFrame.from_dict(distribution_classifications, orient='index', columns=['Correlation', 'Corr_Strength', 'Relation'])
classification_df.index.name = 'Knockdown'
output_path_xlsx = os.path.join(output_dir_classifications, "knockdown_correlation_reasons.xlsx")
classification_df.to_excel(output_path_xlsx)
print (f"Operons: {relation_o}, Pathways: {relation_p}, Regulons: {relation_r}, None: {relation_n}")
print(f"Modneg: {modneg}, noco : {noco}, veweco: {veweco}, vewenegco: {vewenegco}, weco: {weco}, wenegco: {wenegco}, modco: {modco}")
