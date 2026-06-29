# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: PCA Loadings Heatmap (Top Contributing Genes)
# Embedded loading values reconstructed from the screenshot
# For Colab, run: !pip install matplotlib seaborn pandas numpy
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

index = [
    'CTPsyn', 'ACC', 'CPT2', 'bgm', 'Cyt-c', 'Acsl', 'IMPDH2 (CG14806)',
    'ATPsyn-β', 'COX5A', 'Drp1', 'Acon (aconitase)', 'COX4',
    'IMPDH (CG14812)', 'ATPsyn-C', 'Atg8a', 'CRAT', 'Catalase'
]

values = [
    [0.31, -0.05],
    [0.30, -0.02],
    [0.30, -0.09],
    [0.30, -0.08],
    [0.29,  0.12],
    [0.29, -0.15],
    [0.24, -0.24],
    [0.23, -0.05],
    [0.19,  0.21],
    [0.17,  0.24],
    [0.15,  0.37],
    [0.11,  0.45],
    [-0.09, -0.43],
    [-0.14, 0.45],
    [-0.24, 0.22],
    [-0.29, -0.03],
    [-0.29, -0.11]
]

df = pd.DataFrame(values, index=index, columns=['PC1', 'PC2'])

plt.figure(figsize=(11, 10))
sns.heatmap(
    df,
    annot=True,
    fmt='.2f',
    cmap='RdBu_r',
    center=0,
    vmin=-0.43,
    vmax=0.45,
    linewidths=0,
    cbar_kws={'label': 'Loading Value'}
)
plt.title('PCA Loadings Heatmap (Top Contributing Genes)', fontsize=16)
plt.xticks(rotation=0)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig('pca_loadings_heatmap_embedded.png', bbox_inches='tight')
plt.savefig('pca_loadings_heatmap_embedded.pdf', bbox_inches='tight')
plt.show()

df.to_csv('pca_loadings_heatmap_embedded.csv')
print('Saved files: pca_loadings_heatmap_embedded.png, .pdf, and .csv')
