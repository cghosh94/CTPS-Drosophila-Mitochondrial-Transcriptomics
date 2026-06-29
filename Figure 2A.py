# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Distribution of raw p-values / log2FC
# Reconstructed from the screenshot with embedded distributions
# For Colab, run: !pip install matplotlib seaborn pandas numpy
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# p-values: heavy concentration near 0 and 1, broader mid-range, plus edge spikes
pvals = []
pvals.extend(np.random.beta(0.6, 4.5, 1800))
pvals.extend(np.random.uniform(0.05, 0.95, 3300))
pvals.extend(np.random.beta(8, 0.6, 450))
pvals.extend(np.clip(np.random.normal(0.0, 0.01, 220), 0, 1))
pvals.extend(np.clip(np.random.normal(1.0, 0.008, 500), 0, 1))
pvals = np.array(pvals)

# log2FC: narrow center with a long tail
log2fc = np.concatenate([
    np.random.normal(0.0, 0.18, 4200),
    np.random.normal(0.5, 0.35, 700),
    np.random.normal(-0.5, 0.35, 350),
    np.random.normal(1.5, 0.6, 120),
    np.random.normal(-1.5, 0.6, 120),
    np.array([3.2, 3.6, 4.1, -3.3, -4.2, -5.8, 5.1, 6.0])
])

fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.4))

# Left panel: p-values
axes[0].hist(pvals, bins=50, color='#70c7ac', edgecolor='black', alpha=0.95)
axes[0].axvline(0.05, color='#d62728', linestyle='--', linewidth=2.2, label='p=0.05')
axes[0].set_title('Distribution of raw p-values', fontsize=16)
axes[0].set_xlabel('p-value', fontsize=13)
axes[0].set_ylabel('Number of genes', fontsize=13)
axes[0].legend(loc='upper center', frameon=True)
axes[0].set_ylim(0, 530)

# Right panel: log2FC
axes[1].hist(log2fc, bins=40, color='#70c7ac', edgecolor='black', alpha=0.95)
axes[1].set_title('Distribution of log2FC', fontsize=16)
axes[1].set_xlabel('log2 Fold Change', fontsize=13)
axes[1].set_ylabel('Number of genes', fontsize=13)
axes[1].set_xlim(-6.3, 6.5)
axes[1].set_ylim(0, 2100)

for ax in axes:
    ax.grid(True, alpha=0.55)

plt.tight_layout()
plt.savefig('distribution_raw_pvalues_log2fc.png', bbox_inches='tight')
plt.savefig('distribution_raw_pvalues_log2fc.pdf', bbox_inches='tight')
plt.show()

pd.DataFrame({'p_value': pvals}).to_csv('raw_pvalues_embedded.csv', index=False)
pd.DataFrame({'log2FC': log2fc}).to_csv('log2fc_embedded.csv', index=False)
print('Saved files: distribution_raw_pvalues_log2fc.png, .pdf, and .csv files')
