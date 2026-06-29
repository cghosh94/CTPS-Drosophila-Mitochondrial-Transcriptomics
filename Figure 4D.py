# ============================================================
# Correlation heatmap of specified genes in CTPs‑treated samples
# ============================================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 1. Load the count data (adjust the file path if necessary)
# ------------------------------------------------------------------
file_path = 'GSE221707_counts_anno.xls'   # change to full path if needed
try:
    df = pd.read_csv(file_path, sep='\t', index_col=0)
except FileNotFoundError:
    print(f"Error: File '{file_path}' not found. Please check the path.")
    exit()

# The first six columns are the sample counts (as shown in your excerpt)
counts = df.iloc[:, :6]
counts.columns = ['Contr1', 'Contr2', 'Contr3', 'CTPs_Ri1', 'CTPs_Ri2', 'CTPs_Ri3']

# ------------------------------------------------------------------
# 2. Define mapping from common names to possible Drosophila gene symbols
#    (Add synonyms if needed; the script will search case‑insensitively)
# ------------------------------------------------------------------
gene_mapping = {
    'Tom20':   ['Tom20', 'CG7654'],
    'Tom40':   ['Tom40', 'CG12157'],
    'Tom70':   ['Tom70', 'CG6756'],
    'Fis1':    ['Fis1', 'CG3907'],
    'Drp1':    ['Drp1', 'CG32170'],          # Drp1 is present in your data
    'Marf':    ['Marf', 'CG3869'],
    'Miro':    ['Miro', 'CG5410'],
    'Pink1':   ['Pink1', 'CG4523'],
    'Mul1':    ['Mul1', 'CG1134'],
    'Usp30':   ['Usp30', 'CG3016'],
    'Pgam5':   ['Pgam5', 'CG14816'],
    'Tspo':    ['Tspo', 'CG17844'],
    'porin':   ['porin', 'Porin', 'CG6647'],
    'CG7639':  ['CG7639'],
    'CTPS':    ['CTPsyn', 'CG39645']          # CTPS (CTP synthase) – adjust if needed
}

# ------------------------------------------------------------------
# 3. Find which genes are actually present in the data
# ------------------------------------------------------------------
found_genes = {}   # common name -> actual row name in counts
all_index_lower = [idx.lower() for idx in counts.index]

for common, synonyms in gene_mapping.items():
    for syn in synonyms:
        if syn.lower() in all_index_lower:
            actual_name = counts.index[all_index_lower.index(syn.lower())]
            found_genes[common] = actual_name
            break

if not found_genes:
    print("No matching genes found. Check the gene names in your file.")
    print("First 20 gene names in your file:", counts.index[:20].tolist())
    exit()

print(f"Found {len(found_genes)} out of {len(gene_mapping)} requested genes.\n")
for common, actual in found_genes.items():
    print(f"  {common:8} -> {actual}")

# ------------------------------------------------------------------
# 4. Extract expression data for CTPs samples
# ------------------------------------------------------------------
ctp_samples = ['CTPs_Ri1', 'CTPs_Ri2', 'CTPs_Ri3']
available_genes = list(found_genes.values())
ctp_counts = counts.loc[available_genes, ctp_samples].T   # rows = samples, columns = genes

# ------------------------------------------------------------------
# 5. Compute correlation matrix (Pearson)
# ------------------------------------------------------------------
corr_matrix = ctp_counts.corr(method='pearson')
print("\nPearson correlation matrix (CTPs samples):")
print(corr_matrix.round(3))

# ------------------------------------------------------------------
# 6. Plot heatmap
# ------------------------------------------------------------------
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0,
            square=True, linewidths=0.5, cbar_kws={"shrink": 0.8},
            fmt='.2f', annot_kws={'size': 9})
plt.title('Correlation of selected genes\n(CTPs‑treated samples, n=3)', fontsize=14)
plt.tight_layout()
plt.savefig('CTPs_correlation_heatmap.png', dpi=300)
plt.show()

print("\nHeatmap saved as 'CTPs_correlation_heatmap.png'")
