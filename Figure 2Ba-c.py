# ============================================================
# CORRECTED SELF-CONTAINED PYTHON SCRIPT
# Pairwise PCA Scatter Matrix — closer layout/styling match
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Reconstructed sample coordinates so all three pairwise panels match each other
# sample, group, PC1, PC2, PC3
samples = [
    ('C1', 'Control', 4.05,  2.05,  1.50),
    ('C2', 'Control', 3.10,  0.30, -0.80),
    ('C3', 'Control', 1.85, -2.40, -0.60),
    ('R1', 'CTPS_Ri', -3.95,  2.45,  2.02),
    ('R2', 'CTPS_Ri', -4.25, -0.90, -1.45),
    ('R3', 'CTPS_Ri', -0.85, -1.55, -0.68)
]

df = pd.DataFrame(samples, columns=['sample', 'group', 'PC1', 'PC2', 'PC3'])
colors = {'Control': 'blue', 'CTPS_Ri': 'red'}

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
plot_specs = [
    ('PC1', 'PC2', 'PC1 vs PC2', 'PC1 (62.2%)', 'PC2 (19.1%)', axes[0,0]),
    ('PC1', 'PC3', 'PC1 vs PC3', 'PC1 (62.2%)', 'PC3 (9.8%)', axes[0,1]),
    ('PC2', 'PC3', 'PC2 vs PC3', 'PC2 (19.1%)', 'PC3 (9.8%)', axes[1,0])
]

for xcol, ycol, title, xlabel, ylabel, ax in plot_specs:
    for grp in ['Control', 'CTPS_Ri']:
        sub = df[df['group'] == grp]
        ax.scatter(
            sub[xcol], sub[ycol],
            s=55, color=colors[grp], edgecolor='black', linewidth=0.8, alpha=0.95
        )
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.45)

# Match approximate visible axes ranges from screenshot
axes[0,0].set_xlim(-4.7, 4.5)
axes[0,0].set_ylim(-2.6, 2.7)

axes[0,1].set_xlim(-4.7, 4.5)
axes[0,1].set_ylim(-1.6, 2.2)

axes[1,0].set_xlim(-2.6, 2.8)
axes[1,0].set_ylim(-1.6, 2.2)

# Empty panel with centered legend in lower-right quadrant
axes[1,1].axis('off')
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Control', markerfacecolor='blue', markeredgecolor='black', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='CTPS_Ri', markerfacecolor='red', markeredgecolor='black', markersize=8)
]
axes[1,1].legend(handles=legend_elements, loc='center', frameon=True)

fig.suptitle('Pairwise PCA Scatter Matrix', fontsize=18)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('pairwise_pca_scatter_matrix_corrected.png', bbox_inches='tight')
plt.savefig('pairwise_pca_scatter_matrix_corrected.pdf', bbox_inches='tight')
plt.show()

df.to_csv('pairwise_pca_scatter_matrix_corrected.csv', index=False)
print('Saved files: pairwise_pca_scatter_matrix_corrected.png, .pdf, and .csv')
