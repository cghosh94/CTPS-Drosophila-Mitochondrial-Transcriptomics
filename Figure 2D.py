# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: PCA Biplot (Top Contributing Genes)
# Reconstructed from the screenshot with embedded sample scores and gene loadings
# For Colab, run: !pip install matplotlib pandas numpy
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

samples = pd.DataFrame([
    ('Control', 4.0, 2.05),
    ('Control', 3.0, 0.30),
    ('Control', 1.9, -2.40),
    ('CTPS_Ri', -4.0, 2.45),
    ('CTPS_Ri', -4.2, -0.90),
    ('CTPS_Ri', -0.8, -1.55)
], columns=['group', 'PC1', 'PC2'])

loadings = pd.DataFrame([
    ('ATPsyn-C', -0.9,  2.35),
    ('COX4', 0.45,  2.35),
    ('Acon (aconitase)', 0.75, 1.90),
    ('Drp1', 0.95, 1.25),
    ('COX5A', 1.05, 1.05),
    ('Cyt-c', 1.55, 0.65),
    ('ACC', 1.55, -0.10),
    ('CTPsyn', 1.70, -0.18),
    ('ATPsyn-β', 1.30, -0.28),
    ('bgm', 1.25, -0.32),
    ('CPT2', 1.25, -0.42),
    ('Acsl', 1.35, -0.78),
    ('IMPDH2 (CG14806)', 1.35, -1.20),
    ('Atg8a', -1.30, 1.10)
], columns=['gene', 'x', 'y'])

fig, ax = plt.subplots(figsize=(12, 9))

for grp, color in [('Control', 'blue'), ('CTPS_Ri', 'red')]:
    sub = samples[samples['group'] == grp]
    ax.scatter(sub['PC1'], sub['PC2'], s=170, color=color, edgecolor='black', alpha=0.8, label=grp)

for _, row in loadings.iterrows():
    ax.arrow(0, 0, row['x'], row['y'], color='gray', alpha=0.5, width=0.008, head_width=0.07, length_includes_head=True)
    ax.text(row['x']*1.02, row['y']*1.02, row['gene'], fontsize=11)

ax.axhline(0, color='gray', linewidth=0.8)
ax.axvline(0, color='gray', linewidth=0.8)
ax.set_title('PCA Biplot (Top Contributing Genes)', fontsize=16)
ax.set_xlabel('PC1 (62.2%)', fontsize=13)
ax.set_ylabel('PC2 (19.1%)', fontsize=13)
ax.set_xlim(-4.7, 4.5)
ax.set_ylim(-2.7, 2.7)
ax.grid(True, alpha=0.45)
ax.legend(loc='center right', frameon=True)
plt.tight_layout()
plt.savefig('pca_biplot_top_contributing_genes.png', bbox_inches='tight')
plt.savefig('pca_biplot_top_contributing_genes.pdf', bbox_inches='tight')
plt.show()

samples.to_csv('pca_biplot_top_contributing_genes_samples.csv', index=False)
loadings.to_csv('pca_biplot_top_contributing_genes_loadings.csv', index=False)
print('Saved files: pca_biplot_top_contributing_genes.png, .pdf, .csv files')
