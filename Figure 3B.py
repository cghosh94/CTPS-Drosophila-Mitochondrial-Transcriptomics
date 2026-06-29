# ============================================================
# CORRECTED EMBEDDED CODE FOR 3.jpg
# Keeps the screenshot row order exactly and flips the y-axis so
# ACC is on top and Acyl-CoA oxidase is at the bottom, matching
# the reference figure much more closely.
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Row order fixed to match the screenshot from top to bottom
rows = [
    ('ACC', 0.96, 'positive'),
    ('CPT2', 0.95, 'positive'),
    ('Acsl', 0.95, 'positive'),
    ('bgm', 0.95, 'positive'),
    ('Cyt-c', 0.94, 'positive'),
    ('IMPDH2 (CG14806)', 0.83, 'positive'),
    ('HIBCH (CG5044)', 0.82, 'positive'),
    ('ATPsyn-β', 0.73, 'positive'),
    ('COX5A', 0.58, 'positive'),
    ('Drp1', 0.57, 'positive'),
    ('Acon', 0.52, 'positive'),
    ('COX4', 0.32, 'positive'),
    ('SPATA5 (CG5776)', 0.06, 'positive'),
    ('BCS1 (CG4908)', -0.02, 'negative'),
    ('COX18 (CG4942)', -0.10, 'negative'),
    ('CRAT', -0.97, 'negative'),
    ('Catalase', -0.90, 'negative'),
    ('ECI (CG4592)', -0.87, 'negative'),
    ('FUNDC1 (CG5676)', -0.76, 'negative'),
    ('Atg8a', -0.69, 'negative'),
    ('cora (protein 4.1)', -0.66, 'negative'),
    ('MTO1 (CG4610)', -0.63, 'negative'),
    ('ClpP (CG5045)', -0.56, 'negative'),
    ('IDHy (CG5028)', -0.48, 'negative'),
    ('DLST (CG5214)', -0.32, 'negative'),
    ('IMPDH (CG14812)', -0.18, 'negative'),
    ('Ank2', -0.17, 'negative'),
    ('Acyl-CoA oxidase (CG4586)', -0.14, 'negative')
]

df = pd.DataFrame(rows, columns=['Gene', 'Correlation', 'Type'])
colors = {'positive': '#e31a1c', 'negative': '#2a2adf'}

fig, ax = plt.subplots(figsize=(13.0, 7.8))
ax.barh(df['Gene'], df['Correlation'], color=df['Type'].map(colors), edgecolor='white', linewidth=1.0)
ax.invert_yaxis()  # critical correction to match screenshot ordering

ax.set_title('Top 15 Positive and Negative Correlations with CTPS', fontsize=18)
ax.set_xlabel('Pearson correlation with CTPS', fontsize=14)
ax.set_ylabel('Gene', fontsize=14)
ax.set_xlim(-1.08, 1.06)
ax.grid(True, axis='x', alpha=0.45)
ax.grid(False, axis='y')

legend_handles = [
    Patch(color=colors['positive'], label='positive'),
    Patch(color=colors['negative'], label='negative')
]
ax.legend(handles=legend_handles, title='Correlation', loc='lower right', frameon=True)

plt.tight_layout()
plt.savefig('top15_correlations_ctps_corrected.png', bbox_inches='tight')
plt.savefig('top15_correlations_ctps_corrected.pdf', bbox_inches='tight')
plt.show()

df.to_csv('top15_correlations_ctps_corrected.csv', index=False)
print('Saved files: top15_correlations_ctps_corrected.png, .pdf, and .csv')
