# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Differential Expression of All Genes
# No raw file upload needed; embedded gene names and log2FC values used
# For Colab, run: !pip install matplotlib pandas
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Embedded values reconstructed from the plot
# Columns: gene, log2FC, color_group
points = [
    ('CRAT (carnitine acetyltransferase)', 0.70, 'highlight'),
    ('ACAD (CG4860)', 0.66, 'highlight'),
    ('ECI2 (CG4594)', 0.95, 'highlight'),
    ('Serine protease (CG3700)', 2.00, 'highlight'),
    ('CTPsyn', -1.10, 'highlight'),
    ('Cyt-c', -0.30, 'grey'),
    ('PDH-β (CG11876)', -0.66, 'grey'),
    ('carnitine palmitoyltransferase 2', -0.23, 'grey'),
    ('AcCoAS (acetyl-CoA synthetase)', -1.10, 'grey'),
    ('ACC (acetyl-CoA carboxylase)', -0.88, 'grey'),
    ('ATPsyn-α', -0.18, 'grey'),
    ('ECI (CG4592)', 0.43, 'grey'),
    ('IMPDH2 (CG14806)', -0.52, 'grey'),
    ('ATPsyn-β', 0.55, 'grey'),
    ('MTO1 (CG4610)', -0.13, 'grey'),
    ('FATP', 0.24, 'grey'),
    ('ATPSyn-b', 0.13, 'grey'),
    ('MICOS', 0.17, 'grey'),
    ('ClpP', 0.15, 'grey'),
    ('FAM136A', 0.10, 'grey'),
    ('GDT', -0.12, 'grey'),
    ('PRR', 0.16, 'grey'),
    ('COX5A', 0.38, 'grey'),
    ('Drp1', 0.26, 'grey'),
    ('Atg8a', 0.05, 'grey'),
    ('Atg1', 0.20, 'grey'),
    ('ATP synthase subunit', 0.42, 'grey'),
    ('DLST', -0.05, 'grey'),
    ('BCS1', 0.03, 'grey'),
    ('MetaP2', 0.04, 'grey'),
    ('CANT1', 0.02, 'grey'),
    ('COX18', 0.21, 'grey'),
    ('IMPDH', -0.03, 'grey'),
    ('Ank2', 0.01, 'grey')
]

# Add the rest as generic points to preserve the crowded bar layout
more_points = [
    ('g1', -0.10, 'grey'), ('g2', -0.30, 'grey'), ('g3', -0.65, 'grey'), ('g4', -0.22, 'grey'),
    ('g5', -1.10, 'grey'), ('g6', -0.90, 'grey'), ('g7', -1.20, 'grey'), ('g8', 0.43, 'grey'),
    ('g9', -0.52, 'grey'), ('g10', -0.12, 'grey'), ('g11', 0.55, 'grey'), ('g12', 0.24, 'grey'),
    ('g13', 0.14, 'grey'), ('g14', 0.18, 'grey'), ('g15', 0.13, 'grey'), ('g16', 0.11, 'grey'),
    ('g17', -0.06, 'grey'), ('g18', 0.15, 'grey'), ('g19', 0.39, 'grey'), ('g20', -0.16, 'grey'),
    ('g21', -0.09, 'grey'), ('g22', 0.26, 'grey'), ('g23', -0.08, 'grey'), ('g24', -0.05, 'grey'),
    ('g25', 0.52, 'grey'), ('g26', -0.11, 'grey'), ('g27', -0.07, 'grey'), ('g28', -0.03, 'grey'),
    ('g29', -0.05, 'grey'), ('g30', -0.06, 'grey'), ('g31', 0.27, 'grey'), ('g32', 0.13, 'grey'),
    ('g33', -0.12, 'grey'), ('g34', 0.16, 'grey'), ('g35', 0.05, 'grey'), ('g36', 0.17, 'grey'),
    ('g37', 0.42, 'grey'), ('g38', -0.05, 'grey'), ('g39', 0.01, 'grey'), ('g40', 0.20, 'grey'),
    ('g41', 0.14, 'grey'), ('g42', 0.10, 'grey'), ('g43', 0.03, 'grey'), ('g44', 0.27, 'grey'),
    ('g45', -0.18, 'grey'), ('g46', -0.24, 'grey'), ('g47', 0.05, 'grey'), ('g48', 0.06, 'grey'),
    ('g49', -0.03, 'grey'), ('g50', 0.09, 'grey'), ('g51', 0.01, 'grey'), ('g52', 0.02, 'grey'),
    ('g53', 0.10, 'grey'), ('g54', 0.04, 'grey'), ('g55', 0.02, 'grey'), ('g56', 0.07, 'grey'),
    ('g57', -0.04, 'grey'), ('g58', 0.03, 'grey'), ('g59', 0.05, 'grey'), ('g60', 0.42, 'grey'),
    ('g61', -0.02, 'grey'), ('g62', 0.01, 'grey'), ('g63', 0.21, 'grey'), ('g64', 0.02, 'grey'),
    ('g65', 0.00, 'grey'), ('g66', -0.01, 'grey')
]

points = points + more_points

df = pd.DataFrame(points, columns=['gene', 'log2FC', 'group'])

# Sort to mimic the visual order from left to right-ish in the bar chart
# (This is an approximate reconstruction)
df = df.reset_index(drop=True)
df['x'] = range(len(df))

colors = df['group'].map({'highlight': '#e11d21', 'grey': '#8a8a8a'})

plt.figure(figsize=(18, 8))
plt.bar(df['x'], df['log2FC'], color=colors, width=0.8)
plt.axhline(0, color='black', linewidth=1)
plt.title('Differential Expression of All Genes', fontsize=18)
plt.xlabel('Gene', fontsize=12)
plt.ylabel('log2 Fold Change (CTPS_Ri vs Control)', fontsize=12)
plt.xticks(df['x'], df['gene'], rotation=90, fontsize=8)
plt.tight_layout()
plt.savefig('all_genes_differential_expression_embedded.png', bbox_inches='tight')
plt.savefig('all_genes_differential_expression_embedded.pdf', bbox_inches='tight')
plt.show()

df.to_csv('all_genes_differential_expression_embedded.csv', index=False)
print('Saved files: all_genes_differential_expression_embedded.png, .pdf, and .csv')
