# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Top 30 Differentially Expressed Genes heatmap
# Embedded z-score matrix approximated from the uploaded heatmap
# For Colab, run: !pip install matplotlib seaborn pandas numpy
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

samples = ['Contr1', 'Contr2', 'Contr3', 'CTPs_Ri1', 'CTPs_Ri2', 'CTPs_Ri3']

# Embedded matrix visually reconstructed from the figure
z = {
    'CTPsyn':   [1.05, 1.00, 0.75, -0.25, -1.35, -1.20],
    'CG14812':  [-0.95, -0.75, 1.45, -0.30, 1.30, -0.60],
    'CG14806':  [0.75, 0.65, 1.05, -0.50, -0.10, -1.90],
    'ATPsynbeta':[0.10, 1.55, 0.15, 0.55, -1.60, -0.65],
    'ATPsync':  [0.45, -0.20, -1.90, -0.20, 0.40, 1.25],
    'blw':      [1.55, 0.95, 0.10, -0.45, -1.35, -0.60],
    'ATPsynB':  [1.85, 0.75, -0.35, -0.80, -0.35, -0.90],
    'COX4':     [1.95, 0.10, -1.20, -0.85, -0.20, 0.15],
    'COX5A':    [1.90, 0.20, -0.25, -1.00, 0.10, -0.90],
    'COX5B':    [1.75, 0.30, -0.20, -0.35, 0.20, -1.55],
    'COX7A':    [0.30, -0.50, -1.50, -0.20, 1.80, 0.05],
    'Cyt-c1':   [1.95, 0.45, -0.70, -0.20, -0.90, -0.55],
    'Cyt-c-p':  [1.55, 1.05, 0.20, -0.75, -0.95, -0.90],
    'ATPsynD':  [1.85, 0.15, -0.35, -1.10, -0.95, 0.50],
    'ATPsynO':  [1.90, 0.60, -0.65, -0.80, -0.05, -0.90],
    'Acon':     [1.20, 0.30, 0.30, -1.45, -1.15, 0.80],
    'Ald':      [1.10, 1.15, -0.10, -0.15, -1.85, -0.15],
    'CG11876':  [1.30, 1.00, 0.40, -0.10, -1.10, -1.35],
    'CG11899':  [-0.20, -0.10, -1.45, 1.70, 0.65, -0.55],
    'Mdh1':     [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    'Sdha':     [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    'ACC':      [1.25, 0.85, 0.30, 0.05, -1.30, -1.35],
    'Acsl':     [0.85, 0.90, 0.60, 0.25, -1.15, -1.60],
    'bgm':      [1.00, 1.05, 0.30, 0.25, -1.20, -1.40],
    'CRAT':     [-0.90, -1.15, -0.70, 0.85, 1.20, 0.90],
    'CPT2':     [1.05, 0.70, 0.75, 0.15, -1.35, -1.30],
    'Acox57D-d':[0.95, 0.30, 0.95, 0.30, -0.85, -1.70],
    'ATPCL':    [0.90, 0.80, 0.60, 0.30, -0.95, -1.70],
    'AcCoAS':   [1.40, 0.95, 0.20, -0.20, -1.10, -1.25],
    'Drp1':     [0.65, 0.20, -0.25, -0.55, -1.45, 0.95]
}

df = pd.DataFrame(z).T
df.columns = samples

plt.figure(figsize=(8, 10))
sns.heatmap(df, cmap='coolwarm', center=0, vmin=-1.9, vmax=1.9, cbar_kws={'label': 'z-score'})
plt.title('Top 30 Differentially Expressed Genes', fontsize=14)
plt.xticks(rotation=0)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig('top30_differentially_expressed_heatmap_embedded.png', bbox_inches='tight')
plt.savefig('top30_differentially_expressed_heatmap_embedded.pdf', bbox_inches='tight')
plt.show()

df.to_csv('top30_differentially_expressed_heatmap_embedded.csv')
print('Saved files: top30_differentially_expressed_heatmap_embedded.png, .pdf, and .csv')
