# ============================================================
# FIXED EMBEDDED CODE FOR 4.jpg
# Valid linkage matrix, same visual style, no invalid cluster indices
# ============================================================

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

labels = [
    'cora (protein 4.1)', 'Atg8a', 'FUNDC1 (CG5676)', 'ECI (CG4592)', 'CRAT', 'Catalase',
    'Ank2', 'BCS1 (CG4908)', 'DLST (CG5214)', 'IDHy (CG5028)', 'COX18 (CG4942)',
    'ClpP (CG5045)', 'MTO1 (CG4610)', 'Acyl-CoA oxidase (CG4586)', 'IMPDH (CG14812)',
    'COX5A', 'COX4', 'Drp1', 'Acon', 'bgm', 'ACC', 'CPT2', 'Acsl', 'HIBCH (CG5044)',
    'Cyt-c', 'IMPDH2 (CG14806)', 'ATPsyn-β', 'SPATA5 (CG5776)'
]

# Embedded 2D coordinates arranged to reproduce the grouping pattern visually.
# Higher similarity within the upper (blue) and lower (orange) blocks.
X = np.array([
    [0.10, 4.80], [0.35, 4.55], [0.78, 4.30], [1.05, 4.10], [1.55, 3.95], [1.85, 3.80],
    [2.00, 3.55], [2.35, 3.35], [2.70, 3.20], [2.95, 3.00], [2.65, 2.65], [2.80, 2.35],
    [2.95, 2.10], [3.15, 1.95], [3.25, 1.70], [3.45, 0.90], [3.18, 0.78], [2.95, 0.62],
    [2.72, 0.52], [2.48, 0.36], [2.18, 0.22], [1.98, 0.32], [1.84, 0.18], [1.70, 0.04],
    [1.50, 0.12], [1.35, 0.28], [1.20, 0.06], [1.05, -0.10]
])

# Valid linkage derived from these coordinates
Z = linkage(X, method='average', metric='euclidean')

fig, ax = plt.subplots(figsize=(10.8, 8.6))
dendrogram(
    Z,
    labels=labels,
    orientation='left',
    leaf_font_size=9,
    color_threshold=0.7 * Z[:, 2].max(),
    above_threshold_color='#8da0cb',
    ax=ax
)

ax.set_title('Hierarchical Clustering of Mitochondrial Genes (based on correlation)', fontsize=14)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=8)
plt.tight_layout()
plt.savefig('hierarchical_clustering_mito_genes_fixed.png', bbox_inches='tight')
plt.savefig('hierarchical_clustering_mito_genes_fixed.pdf', bbox_inches='tight')
plt.show()

pd.DataFrame(X, columns=['x', 'y'], index=labels).to_csv('hierarchical_clustering_mito_genes_embedded_points.csv')
pd.DataFrame(Z, columns=['idx1', 'idx2', 'distance', 'count']).to_csv('hierarchical_clustering_mito_genes_valid_linkage.csv', index=False)
print('Saved files: hierarchical_clustering_mito_genes_fixed.png, .pdf, embedded_points.csv, valid_linkage.csv')
