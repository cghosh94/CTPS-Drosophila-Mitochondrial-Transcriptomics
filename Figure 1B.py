# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Organelle category stacked bars
# Reconstructed from the screenshot with embedded counts/percentages
# For Colab, run: !pip install matplotlib pandas numpy
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

organelles = ['Mitochondria', 'ER', 'Golgi', 'Nucleus', 'Lipid droplet', 'Peroxisome']
counts = [1645, 654, 331, 2127, 49, 77]
up = [12.5, 10.0, 12.0, 10.0, 8.0, 10.5]
down = [3.5, 4.5, 1.0, 1.0, 12.0, 0.8]
unchanged = [100 - u - d for u, d in zip(up, down)]

df = pd.DataFrame({'Organelle': organelles, 'Up': up, 'Down': down, 'Unchanged': unchanged, 'n': counts})

fig, ax = plt.subplots(figsize=(14.6, 8))
bar_x = range(len(organelles))
width = 0.8

ax.bar(bar_x, df['Up'], color='#2ca02c', width=width, label='Up-regulated (p<0.05)')
ax.bar(bar_x, df['Down'], bottom=df['Up'], color='#d62728', width=width, label='Down-regulated (p<0.05)')
ax.bar(bar_x, df['Unchanged'], bottom=df['Up'] + df['Down'], color='lightgray', width=width, label='Unchanged')

for i, n in enumerate(counts):
    ax.text(i, 102, f'n={n}', ha='center', va='bottom', fontsize=13)

ax.set_xticks(list(bar_x))
ax.set_xticklabels(organelles, rotation=45, ha='right', fontsize=12)
ax.set_ylabel('Percentage of genes (%)', fontsize=14)
ax.set_xlabel('Organelle', fontsize=14)
ax.set_ylim(0, 112)
ax.set_title('Mitochondrial genes show the highest proportion of significantly changing genes (p<0.05, any FC)', fontsize=16)
ax.legend(loc='center right', frameon=True, fontsize=12)
ax.grid(True, axis='y', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig('organelle_stacked_bar_embedded.png', bbox_inches='tight')
plt.savefig('organelle_stacked_bar_embedded.pdf', bbox_inches='tight')
plt.show()

df.to_csv('organelle_stacked_bar_embedded.csv', index=False)
print('Saved files: organelle_stacked_bar_embedded.png, .pdf, and .csv')
