# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Scree Plot
# Embedded variance values reconstructed from the screenshot
# For Colab, run: !pip install matplotlib pandas numpy
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

pcs = [1, 2, 3, 4, 5, 6]
individual = [62.2, 19.1, 9.8, 7.0, 1.9, 0.0]
cumulative = [62.2, 81.3, 91.1, 98.1, 100.0, 100.0]

df = pd.DataFrame({'PC': pcs, 'Individual': individual, 'Cumulative': cumulative})

fig, ax = plt.subplots(figsize=(10, 6.5))
ax.bar(df['PC'], df['Individual'], color='#7ea7c9', alpha=0.9, label='Individual')
ax.plot(df['PC'], df['Cumulative'], color='red', marker='o', linewidth=2.5, markersize=8, label='Cumulative')

for x, y in zip(df['PC'], df['Individual']):
    ax.text(x, y + 1.2, f'{y:.1f}%', ha='center', va='bottom', fontsize=11)

ax.set_title('Scree Plot', fontsize=16)
ax.set_xlabel('Principal Component', fontsize=13)
ax.set_ylabel('Variance Explained (%)', fontsize=13)
ax.set_xticks(df['PC'])
ax.set_ylim(0, 105)
ax.legend(loc='upper left', frameon=True)
ax.grid(True, axis='both', alpha=0.5)
plt.tight_layout()
plt.savefig('scree_plot_embedded.png', bbox_inches='tight')
plt.savefig('scree_plot_embedded.pdf', bbox_inches='tight')
plt.show()

df.to_csv('scree_plot_embedded.csv', index=False)
print('Saved files: scree_plot_embedded.png, .pdf, and .csv')
