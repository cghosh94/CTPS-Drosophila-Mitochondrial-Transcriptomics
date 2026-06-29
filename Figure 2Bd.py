# ============================================================
# CORRECTED SELF-CONTAINED PYTHON SCRIPT
# 3D PCA Plot (PC1, PC2, PC3) — closer to screenshot
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.lines import Line2D

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Coordinates adjusted to match the visible screenshot labels/positions more closely
samples = pd.DataFrame([
    ('Contr1',   'Control',  4.10,  2.25,  1.35),
    ('Contr2',   'Control',  2.15, -0.85, -0.95),
    ('Contr3',   'Control',  0.95, -1.85, -1.45),
    ('CTPs_Ri1', 'CTPS_Ri', -1.25, -0.95, -1.20),
    ('CTPs_Ri2', 'CTPS_Ri', -3.95, -2.25,  2.05),
    ('CTPs_Ri3', 'CTPS_Ri', -1.10, -0.35, -0.05)
], columns=['sample', 'group', 'PC1', 'PC2', 'PC3'])

colors = {'Control': 'blue', 'CTPS_Ri': 'red'}

fig = plt.figure(figsize=(10.5, 8))
ax = fig.add_subplot(111, projection='3d')

for grp in ['Control', 'CTPS_Ri']:
    sub = samples[samples['group'] == grp]
    ax.scatter(
        sub['PC1'], sub['PC2'], sub['PC3'],
        s=95, color=colors[grp], edgecolor='black', linewidth=1, alpha=0.8
    )
    for _, row in sub.iterrows():
        ax.text(row['PC1'], row['PC2'], row['PC3'], row['sample'], fontsize=10)

ax.set_title('3D PCA Plot (PC1, PC2, PC3)', fontsize=18)
ax.set_xlabel('PC1 (62.2%)', fontsize=12, labelpad=10)
ax.set_ylabel('PC2 (19.1%)', fontsize=12, labelpad=10)
ax.set_zlabel('PC3 (9.8%)', fontsize=12, labelpad=8)

# Match general camera angle / box feel from screenshot
ax.view_init(elev=31, azim=-61)
ax.set_xlim(-4.5, 4.5)
ax.set_ylim(-2.7, 2.7)
ax.set_zlim(-1.6, 2.2)

legend_handles = [
    Line2D([0], [0], marker='o', color='w', label='Control', markerfacecolor='blue', markeredgecolor='white', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='CTPS_Ri', markerfacecolor='red', markeredgecolor='white', markersize=8)
]
ax.legend(handles=legend_handles, loc='upper left', frameon=True)

plt.tight_layout()
plt.savefig('3d_pca_plot_corrected.png', bbox_inches='tight')
plt.savefig('3d_pca_plot_corrected.pdf', bbox_inches='tight')
plt.show()

samples.to_csv('3d_pca_plot_corrected.csv', index=False)
print('Saved files: 3d_pca_plot_corrected.png, .pdf, and .csv')
