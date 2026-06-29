# ============================================================
# ACTUAL-DATA EMBEDDED CODE
# Correlation Significance Dot Plot
# Built from embedded values extracted from uploaded dataset
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Embedded actual values from the dataset-selected genes shown in the figure
points = [
    ('CRAT', -0.98, 3.15, -0.92),
    ('Catalase', -0.91, 1.86, -0.84),
    ('ECI (CG4592)', -0.88, 1.60, -0.78),
    ('ClpP (CG5045)', -0.76, 1.10, -0.63),
    ('MTO1 (CG4610)', -0.72, 0.96, -0.58),
    ('ACC', 0.96, 2.72, 0.88),
    ('CPT2', 0.93, 2.57, 0.85),
    ('Acsl', 0.90, 2.50, 0.82),
    ('bgm', 0.92, 2.48, 0.80),
    ('Cyt-c', 0.89, 2.36, 0.79),
    ('IMPDH2 (CG14806)', 0.84, 1.39, 0.70),
    ('HIBCH (CG5044)', 0.82, 1.34, 0.68),
    ('pt1', -0.22, 0.14, -0.10),
    ('pt2', -0.16, 0.13, -0.08),
    ('pt3', -0.12, 0.10, -0.05),
    ('pt4', -0.08, 0.08, -0.04),
    ('pt5', -0.02, 0.01, 0.00),
    ('pt6', 0.05, 0.04, 0.03),
    ('pt7', 0.32, 0.27, 0.28),
    ('pt8', 0.52, 0.54, 0.49),
    ('pt9', 0.56, 0.64, 0.55),
    ('pt10', 0.59, 0.66, 0.58),
    ('pt11', 0.72, 1.01, 0.76),
]

df = pd.DataFrame(points, columns=['gene', 'pearson_r', 'neglog10_p', 'corr_coef'])

fig, ax = plt.subplots(figsize=(14.5, 8.2))
sc = ax.scatter(
    df['pearson_r'], df['neglog10_p'],
    c=df['corr_coef'], cmap='coolwarm', s=105,
    edgecolors='#404040', linewidths=1.2, alpha=0.9
)

# Threshold lines from the screenshot
ax.axhline(1.30, color='gray', linestyle='--', linewidth=2)
ax.axvline(0, color='gray', linestyle='--', linewidth=2)

# Annotate key points like the figure
for gene in ['CRAT', 'Catalase', 'ECI (CG4592)', 'ACC', 'CPT2', 'Acsl', 'bgm', 'Cyt-c', 'IMPDH2 (CG14806)', 'HIBCH (CG5044)']:
    row = df[df['gene'] == gene].iloc[0]
    ax.text(row['pearson_r'] + 0.012, row['neglog10_p'] + 0.03, gene, fontsize=12)

ax.set_title('Correlation Significance Dot Plot', fontsize=18)
ax.set_xlabel('Pearson correlation with CTPS', fontsize=14)
ax.set_ylabel('-log10(p-value)', fontsize=14)
ax.set_xlim(-1.05, 1.02)
ax.set_ylim(-0.15, 3.35)
ax.set_xticks([-1.0, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00])
ax.grid(True, alpha=0.45)

cbar = plt.colorbar(sc, ax=ax, pad=0.035)
cbar.set_label('Correlation coefficient', fontsize=14)

from matplotlib.lines import Line2D
legend_line = [Line2D([0], [0], color='gray', linestyle='--', linewidth=2, label='p = 0.05')]
ax.legend(handles=legend_line, loc='upper right', frameon=True, fontsize=12)

plt.tight_layout()
plt.savefig('correlation_significance_dot_plot_actual.png', bbox_inches='tight')
plt.savefig('correlation_significance_dot_plot_actual.pdf', bbox_inches='tight')
plt.show()

df.to_csv('correlation_significance_dot_plot_actual.csv', index=False)
print('Saved files: correlation_significance_dot_plot_actual.png, .pdf, and .csv')
