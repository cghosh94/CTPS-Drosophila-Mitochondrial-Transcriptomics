# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Volcano plot: Mitochondrial genes (Control vs CTPs)
# Reconstructed from the screenshot with embedded points and highlighted genes
# For Colab, run: !pip install matplotlib pandas numpy
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

np.random.seed(7)

# Background mitochondrial genes: dense center + spread outliers to mimic the screenshot
bg = []
for i in range(150):
    x = np.random.normal(0.0, 0.35)
    y = abs(np.random.normal(0.8, 0.55))
    bg.append((f'bg_{i}', x, min(y, 4.2), 'Other mitochondrial genes'))
for i in range(70):
    x = np.random.normal(0.25, 0.65)
    y = abs(np.random.normal(1.4, 0.75))
    bg.append((f'bg2_{i}', x, min(y, 4.2), 'Other mitochondrial genes'))
for i, (x, y) in enumerate([(-3.9, 0.62), (-3.2, 0.92), (-2.8, 0.62), (-1.9, 0.52), (-1.6, 1.23), (1.55, 4.18), (2.18, 4.12), (1.65, 2.82), (2.72, 1.28)]):
    bg.append((f'out_{i}', x, y, 'Other mitochondrial genes'))

# Highlighted outer membrane genes from the figure
highlight = [
    ('Fis1', 0.26, 1.95, 'Outer membrane genes'),
    ('Tom20', 0.14, 1.38, 'Outer membrane genes'),
    ('Tom40', -0.08, 0.66, 'Outer membrane genes'),
    ('Pgam5', 0.20, 0.72, 'Outer membrane genes'),
    ('Tom70', 0.03, 0.18, 'Outer membrane genes'),
    ('Pink1', 0.05, 0.12, 'Outer membrane genes'),
    ('Mul1', -0.01, 0.01, 'Outer membrane genes'),
    ('Miro', -0.02, 0.05, 'Outer membrane genes')
]

all_points = bg + highlight

df = pd.DataFrame(all_points, columns=['gene', 'log2FC', 'neglog10_pvalue', 'group'])

fig, ax = plt.subplots(figsize=(10, 8))

other = df[df['group'] == 'Other mitochondrial genes']
outm = df[df['group'] == 'Outer membrane genes']

ax.scatter(other['log2FC'], other['neglog10_pvalue'], color='lightgray', s=30, alpha=0.8, label='Other mitochondrial genes')
ax.scatter(outm['log2FC'], outm['neglog10_pvalue'], color='red', edgecolor='black', s=85, label='Outer membrane genes', zorder=3)

# Threshold lines as shown in the panel
ax.axhline(-np.log10(0.05), color='#4f67ff', linestyle='--', linewidth=1.6, alpha=0.7, label='p = 0.05')
ax.axvline(-1, color='#4f67ff', linestyle='--', linewidth=1.6, alpha=0.7)
ax.axvline(1, color='#4f67ff', linestyle='--', linewidth=1.6, alpha=0.7)

for _, row in outm.iterrows():
    ax.text(row['log2FC'] + 0.03, row['neglog10_pvalue'] + 0.03, row['gene'], fontsize=11)

ax.set_title('Volcano plot: Mitochondrial genes (Control vs CTPs)', fontsize=16)
ax.set_xlabel('log$_2$ fold change (CTPs / Control)')
ax.set_ylabel('-log$_{10}$(p-value)')
ax.set_xlim(-4.2, 2.9)
ax.set_ylim(-0.2, 4.4)
ax.legend(loc='upper left', frameon=True)
plt.tight_layout()
plt.savefig('volcano_mito_outer_membrane_control_vs_ctps.png', bbox_inches='tight')
plt.savefig('volcano_mito_outer_membrane_control_vs_ctps.pdf', bbox_inches='tight')
plt.show()

df.to_csv('volcano_mito_outer_membrane_control_vs_ctps.csv', index=False)
print('Saved files: volcano_mito_outer_membrane_control_vs_ctps.png, .pdf, and .csv')
