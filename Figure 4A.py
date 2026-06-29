# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Mitochondrial outer membrane genes: Control vs CTPs
# Reconstructed from the screenshot with embedded means, errors, and p-values
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

records = [
    ('Tom20', 'Control', 5950, 100, 'p=0.041'), ('Tom20', 'CTPs', 6500, 110, 'p=0.041'),
    ('Tom40', 'Control', 9200, 130, 'p=0.220'), ('Tom40', 'CTPs', 8700, 120, 'p=0.220'),
    ('Tom70', 'Control', 2350, 60, 'p=0.654'), ('Tom70', 'CTPs', 2400, 55, 'p=0.654'),
    ('Fis1', 'Control', 650, 25, 'p=0.012'), ('Fis1', 'CTPs', 820, 30, 'p=0.012'),
    ('Drp1', 'Control', 1800, 70, 'p=0.236'), ('Drp1', 'CTPs', 1700, 65, 'p=0.236'),
    ('Marf', 'Control', 2900, 55, 'p=0.958'), ('Marf', 'CTPs', 3000, 50, 'p=0.958'),
    ('Miro', 'Control', 1050, 40, 'p=0.932'), ('Miro', 'CTPs', 980, 35, 'p=0.932'),
    ('Pink1', 'Control', 1250, 40, 'p=0.875'), ('Pink1', 'CTPs', 1220, 35, 'p=0.875'),
    ('Mul1', 'Control', 280, 12, 'p=0.657'), ('Mul1', 'CTPs', 240, 10, 'p=0.657'),
    ('Usp30', 'Control', 650, 20, 'p=0.246'), ('Usp30', 'CTPs', 720, 22, 'p=0.246'),
    ('Pgam5', 'Control', 780, 25, 'p=0.189'), ('Pgam5', 'CTPs', 980, 28, 'p=0.189'),
    ('Tspo', 'Control', 650, 18, 'p=0.083'), ('Tspo', 'CTPs', 1200, 30, 'p=0.083'),
    ('porin', 'Control', 43500, 900, 'p=0.078'), ('porin', 'CTPs', 46000, 800, 'p=0.078'),
    ('CG7639', 'Control', 1400, 45, 'p=0.181'), ('CG7639', 'CTPs', 1550, 50, 'p=0.181')
]

df = pd.DataFrame(records, columns=['Gene', 'group', 'expression', 'err', 'p'])
order = ['Control', 'CTPs']
palette = {'Control': '#66c2a5', 'CTPs': '#fc8d62'}
genes = ['Tom20', 'Tom40', 'Tom70', 'Fis1', 'Drp1', 'Marf', 'Miro', 'Pink1', 'Mul1', 'Usp30', 'Pgam5', 'Tspo', 'porin', 'CG7639']

fig, ax = plt.subplots(figsize=(19, 8))
bar_width = 0.35
x = range(len(genes))

for i, group in enumerate(order):
    sub = df[df['group'] == group].set_index('Gene').loc[genes].reset_index()
    xpos = [j + (i - 0.5) * bar_width for j in x]
    ax.bar(xpos, sub['expression'], width=bar_width, color=palette[group], label=group, edgecolor='none')
    ax.errorbar(xpos, sub['expression'], yerr=sub['err'], fmt='none', ecolor='black', elinewidth=2, capsize=6, capthick=1.2)

for i, gene in enumerate(genes):
    p = df[df['Gene'] == gene]['p'].iloc[0]
    y = df[df['Gene'] == gene]['expression'].max() + df[df['Gene'] == gene]['err'].max() + (250 if gene != 'porin' else 700)
    if gene == 'porin':
        y = 46800
    ax.text(i, y, p, ha='center', va='bottom', fontsize=12)

ax.set_xticks(list(x))
ax.set_xticklabels(genes, rotation=45, ha='right')
ax.set_ylabel('Normalized expression')
ax.set_title('Mitochondrial outer membrane genes: Control vs CTPs')
ax.legend(loc='upper right', frameon=True)
ax.set_ylim(0, 48000)
plt.tight_layout()
plt.savefig('mitochondrial_outer_membrane_control_vs_ctps.png', bbox_inches='tight')
plt.savefig('mitochondrial_outer_membrane_control_vs_ctps.pdf', bbox_inches='tight')
plt.show()

df.to_csv('mitochondrial_outer_membrane_control_vs_ctps.csv', index=False)
print('Saved files: mitochondrial_outer_membrane_control_vs_ctps.png, .pdf, and .csv')
