# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Expression Changes by Functional Category
# 'Other' category moved to the last position
# For Colab, run: !pip install matplotlib seaborn pandas numpy
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style='whitegrid', context='talk')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

records = []

def add(cat, vals):
    for v in vals:
        records.append((cat, v))

add('Fatty acid metabolism', [-1.24, -1.10, -0.98, -0.80, -0.72, -0.60, -0.35, -0.25, 0.72])
add('Other', [-1.10, -0.52, -0.25, -0.18, -0.05, 0.05, 0.12, 0.18, 0.42, 0.66, 0.95, 2.18])
add('Stress/Autophagy', [0.22, 0.28, 0.32, 0.38, 0.46, 0.50, 0.55, 2.01])
add('ATP synthesis', [-0.30, -0.18, -0.14, -0.12, -0.10, -0.08, -0.05, 0.06])
add('TCA cycle', [-0.68, -0.22, -0.18, -0.12, -0.10, -0.05, 0.00, 0.12, 0.26])
add('Mitochondrial dynamics', [-0.12, -0.10, -0.08, -0.06, -0.02, 0.00, 0.06, 0.12, 0.20, 0.42])

df = pd.DataFrame(records, columns=['Category', 'log2FC'])

order = [
    'Fatty acid metabolism',
    'Stress/Autophagy',
    'ATP synthesis',
    'TCA cycle',
    'Mitochondrial dynamics',
    'Other'
]

palette = ['#8ecfc9', '#beb8d8', '#f39c97', '#89b0d3', '#f1b56b', '#f6f5a6']

fig, ax = plt.subplots(figsize=(15, 8))
sns.boxplot(
    data=df,
    x='Category',
    y='log2FC',
    order=order,
    palette=palette,
    width=0.75,
    fliersize=9,
    linewidth=1.5,
    ax=ax
)

ax.axhline(0, color='black', linestyle='--', linewidth=2)
ax.set_title('Expression Changes by Functional Category', fontsize=18)
ax.set_xlabel('Category', fontsize=14)
ax.set_ylabel('log2 Fold Change', fontsize=14)
ax.set_ylim(-1.35, 2.30)
ax.tick_params(axis='x', rotation=0)

plt.tight_layout()
plt.savefig('functional_category_boxplot_other_last.png', bbox_inches='tight')
plt.savefig('functional_category_boxplot_other_last.pdf', bbox_inches='tight')
plt.show()

df.to_csv('functional_category_boxplot_other_last.csv', index=False)
print('Saved files: functional_category_boxplot_other_last.png, .pdf, and .csv')
