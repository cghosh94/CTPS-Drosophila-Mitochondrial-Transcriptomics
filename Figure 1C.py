# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Heatmap of top 5 most significant genes per organelle
# Reconstructed from the screenshot with embedded normalized values
# For Colab, run: !pip install matplotlib seaborn pandas numpy scipy
# ============================================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

samples = ['Contr3', 'Contr1', 'Contr2', 'CTPs_Ri1', 'CTPs_Ri2', 'CTPs_Ri3']

rows = [
    ('ProRS-m', 'Mitochondria', [0.05, 0.05, 0.04, 0.92, 0.90, 0.88]),
    ('CG12264', 'Mitochondria', [0.04, 0.12, 0.08, 0.70, 0.75, 0.90]),
    ('CG30022', 'Mitochondria', [0.04, 0.10, 0.12, 0.72, 0.95, 0.98]),
    ('Gnf5', 'Mitochondria', [0.10, 0.10, 0.04, 0.74, 0.72, 0.96]),
    ('CG32581', 'Mitochondria', [0.04, 0.10, 0.04, 0.70, 0.65, 0.95]),

    ('CG1717', 'ER', [0.96, 0.96, 0.96, 0.04, 0.04, 0.04]),
    ('CG31021', 'ER', [0.06, 0.05, 0.05, 0.93, 0.94, 0.90]),
    ('Cyp6d2', 'ER', [0.08, 0.07, 0.05, 0.85, 0.95, 0.86]),
    ('CG15661', 'ER', [0.05, 0.04, 0.10, 0.82, 0.88, 0.95]),
    ('BI-1', 'ER', [0.04, 0.05, 0.10, 0.78, 0.95, 0.95]),

    ('FucTB', 'Golgi', [0.10, 0.12, 0.06, 0.90, 0.72, 0.95]),
    ('obe', 'Golgi', [0.06, 0.10, 0.07, 0.72, 0.72, 0.96]),
    ('Fur1', 'Golgi', [0.12, 0.10, 0.06, 0.80, 0.65, 0.95]),
    ('Ect3', 'Golgi', [0.05, 0.10, 0.14, 0.72, 0.75, 0.95]),
    ('CG31019', 'Golgi', [0.12, 0.05, 0.05, 0.72, 0.95, 0.90]),

    ('HMGZ', 'Nucleus', [0.96, 0.96, 0.96, 0.04, 0.04, 0.04]),
    ('murs30', 'Nucleus', [0.06, 0.05, 0.04, 0.72, 0.80, 0.98]),
    ('Prosbeta4', 'Nucleus', [0.06, 0.12, 0.04, 0.78, 0.75, 0.96]),
    ('Obp', 'Nucleus', [0.10, 0.05, 0.04, 0.72, 0.70, 0.96]),
    ('tyr', 'Nucleus', [0.08, 0.10, 0.04, 0.78, 0.80, 0.96]),

    ('dob', 'Lipid droplet', [0.10, 0.15, 0.04, 0.95, 0.65, 0.65]),
    ('CG5112', 'Lipid droplet', [0.12, 0.12, 0.04, 0.65, 0.90, 0.95]),
    ('Eno', 'Peroxisome', [0.55, 0.96, 0.95, 0.35, 0.04, 0.12]),
    ('Gapdh1', 'Peroxisome', [0.72, 0.95, 0.96, 0.30, 0.04, 0.10]),
    ('Pgi', 'Peroxisome', [0.08, 0.10, 0.04, 0.65, 0.72, 0.96]),
    ('CG12428', 'Peroxisome', [0.15, 0.35, 0.04, 0.62, 0.75, 0.95]),
    ('CG7970', 'Peroxisome', [0.05, 0.30, 0.05, 0.85, 0.95, 0.70]),
    ('CG11737', 'Peroxisome', [0.30, 0.25, 0.20, 0.96, 0.96, 0.45]),
    ('Pex19', 'Peroxisome', [0.60, 0.95, 0.65, 0.05, 0.04, 0.08]),
    ('CG13827', 'Peroxisome', [0.35, 0.05, 0.28, 0.10, 0.96, 0.96]),
]

genes = [r[0] for r in rows]
organelles = [r[1] for r in rows]
expr = np.array([r[2] for r in rows])
heat_df = pd.DataFrame(expr, index=genes, columns=samples)

org_palette = {
    'Mitochondria': '#66c2a5',
    'ER': '#fc8d62',
    'Golgi': '#8da0cb',
    'Nucleus': '#e78ac3',
    'Lipid droplet': '#a6d854',
    'Peroxisome': '#ffd92f'
}
row_colors = pd.Series(organelles, index=genes).map(org_palette)

cg = sns.clustermap(
    heat_df,
    cmap='coolwarm',
    vmin=0,
    vmax=1,
    row_cluster=False,
    col_cluster=True,
    row_colors=row_colors,
    linewidths=0,
    figsize=(8, 10),
    cbar_pos=(0.02, 0.86, 0.12, 0.07),
    cbar_kws={'orientation': 'horizontal'}
)

cg.ax_heatmap.set_title('Top 5 most significant genes per organelle\n(CTPsyn knockdown vs control)', pad=150, fontsize=14)
cg.ax_heatmap.set_xlabel('')
cg.ax_heatmap.set_ylabel('gene_id')
cg.ax_heatmap.tick_params(axis='x', labelrotation=0, labelsize=8)
cg.ax_heatmap.tick_params(axis='y', labelsize=6)

# Label row color strip
cg.ax_row_colors.set_ylabel('Organelle', rotation=90, fontsize=8)

# Legend for organelles
for label, color in org_palette.items():
    cg.ax_heatmap.bar(0, 0, color=color, label=label, linewidth=0)
legend = cg.ax_heatmap.legend(title='Organelle', loc='upper left', bbox_to_anchor=(1.10, 1.0), frameon=True, fontsize=7, title_fontsize=8)

plt.savefig('top5_significant_genes_heatmap_embedded.png', bbox_inches='tight')
plt.savefig('top5_significant_genes_heatmap_embedded.pdf', bbox_inches='tight')
plt.show()

heat_df.insert(0, 'Organelle', organelles)
heat_df.to_csv('top5_significant_genes_heatmap_embedded.csv')
print('Saved files: top5_significant_genes_heatmap_embedded.png, .pdf, and .csv')
