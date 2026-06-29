import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

# Data
ctpsyn = [3777, 3780, 3400, 2211, 1403, 1465]
ras    = [5581, 5217, 4533, 7444, 6344, 5995]

data = pd.DataFrame({
    'Sample': ['control1','control2','control3','CTPS_Ri1','CTPS_Ri2','CTPS_Ri3'],
    'CTPsyn': ctpsyn,
    'ras': ras
})

data_melt = data.melt(id_vars='Sample', var_name='Gene', value_name='Expression')
data_melt['Condition'] = data_melt['Sample'].apply(lambda x: 'Control' if 'control' in x else 'CTPS_Ri')

# Compute means and SEM
stats_df = data_melt.groupby(['Gene', 'Condition'])['Expression'].agg(['mean', 'sem']).reset_index()

condition_order = ['Control', 'CTPS_Ri']
gene_order = ['CTPsyn', 'ras']

# t-tests
p_vals = {}
for gene in gene_order:
    ctrl = data_melt[(data_melt['Gene'] == gene) & (data_melt['Condition'] == 'Control')]['Expression']
    ri   = data_melt[(data_melt['Gene'] == gene) & (data_melt['Condition'] == 'CTPS_Ri')]['Expression']
    _, p = stats.ttest_ind(ctrl, ri)
    p_vals[gene] = p

def star(p):
    if p < 0.0001: return '****'
    if p < 0.001:  return '***'
    if p < 0.01:   return '**'
    if p < 0.05:   return '*'
    return 'ns'

# Plot
plt.figure(figsize=(6,6))
sns.set_style("whitegrid")
sns.set_palette("pastel")

ax = sns.barplot(data=stats_df, x='Gene', y='mean', hue='Condition',
                 hue_order=condition_order, order=gene_order,
                 errorbar=('ci', 95), capsize=0.1, errwidth=1.5, alpha=0.7)

# Overlay individual points with slight jitter to see all three
sns.stripplot(data=data_melt, x='Gene', y='Expression', hue='Condition',
              hue_order=condition_order, order=gene_order,
              dodge=True, jitter=0.15, size=7, edgecolor='black', linewidth=1,
              alpha=0.9, palette=['#FFB3BA', '#BAFFC9'])

# Significance: place horizontal line and star just above the taller bar for each gene
for i, gene in enumerate(gene_order):
    # Get the mean values for Control and CTPS_Ri for this gene
    means = stats_df[stats_df['Gene'] == gene]['mean'].values
    max_mean = max(means)
    sig_y = max_mean + max_mean * 0.05  # 5% above the tallest bar
    
    ax.plot([i-0.2, i+0.2], [sig_y, sig_y], color='black', lw=1.5)
    ax.text(i, sig_y + max_mean*0.01, star(p_vals[gene]), ha='center', va='bottom', fontsize=12, fontweight='bold')

# Legend: upper left corner, remove duplicate entries
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], title='Condition', loc='upper left')

plt.ylabel('Expression level', fontsize=12)
plt.title("CTPsyn and ras expression in Drosophila\nControl vs CTPS_Ri knockdown", fontsize=13)

# Adjust y-limit to leave a little room above the highest star
y_max_plot = max([stats_df[stats_df['Gene']==g]['mean'].max() for g in gene_order])
plt.ylim(0, y_max_plot * 1.15)

plt.tight_layout()
plt.show()
