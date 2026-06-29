# ============================================================
# SELF-CONTAINED PYTHON SCRIPT
# Figure: Organelle-wise volcano-style scatter panels
# Reconstructed from the screenshot with embedded synthetic coordinates
# For Colab, run: !pip install matplotlib pandas numpy
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.random.seed(10)
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 120
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# helper to make panel-like clouds
def make_panel(n, x_mu=0.2, x_sd=0.5, y_scale=1.0, x_clip=(-4,4), y_clip=(0,10), outliers=None):
    x = np.random.normal(x_mu, x_sd, n)
    y = np.random.gamma(shape=1.8, scale=0.55*y_scale, size=n)
    x = np.clip(x, x_clip[0], x_clip[1])
    y = np.clip(y, y_clip[0], y_clip[1])
    if outliers:
        xo, yo = zip(*outliers)
        x = np.concatenate([x, np.array(xo)])
        y = np.concatenate([y, np.array(yo)])
    return pd.DataFrame({'log2FC': x, 'neglog10p': y})

panels = {
    'Mitochondria (n=1645)': make_panel(700, x_mu=0.15, x_sd=0.55, y_scale=1.15, x_clip=(-3.8,3.8), y_clip=(0,4.0), outliers=[(-3.0,1.6),(3.8,3.4),(2.9,2.3),(-2.3,1.9)]),
    'ER (n=654)': make_panel(180, x_mu=0.1, x_sd=0.75, y_scale=0.95, x_clip=(-6,6), y_clip=(0,10), outliers=[(-5.7,3.6),(-5.3,2.5),(-1.0,10.0),(1.1,10.0),(5.8,5.8),(3.1,5.2)]),
    'Golgi (n=331)': make_panel(150, x_mu=0.35, x_sd=0.45, y_scale=0.95, x_clip=(-1.8,4.0), y_clip=(0,3.6), outliers=[(4.0,0.45),(3.2,0.95),(1.8,3.1),(0.5,3.2)]),
    'Nucleus (n=2127)': make_panel(500, x_mu=0.35, x_sd=0.75, y_scale=0.9, x_clip=(-3,6), y_clip=(0,10), outliers=[(-1.0,10.0),(5.9,3.7),(4.0,3.3)]),
    'Lipid droplet (n=49)': make_panel(35, x_mu=-0.1, x_sd=0.8, y_scale=0.7, x_clip=(-3.6,1.6), y_clip=(0,2.6), outliers=[(-3.5,0.7),(-2.6,1.45),(1.2,2.55),(0.9,2.4)]),
    'Peroxisome (n=77)': make_panel(50, x_mu=0.15, x_sd=0.55, y_scale=0.7, x_clip=(-1.5,2.2), y_clip=(0,3.0), outliers=[(2.2,2.98),(1.9,0.82),(-0.6,1.9)])
}

fig, axes = plt.subplots(2, 3, figsize=(14.8, 9.2))
axes = axes.flatten()
threshold = -np.log10(0.05)
color = '#66c2a4'

for ax, (title, df) in zip(axes, panels.items()):
    ax.scatter(df['log2FC'], df['neglog10p'], s=20, color=color, alpha=0.6, edgecolors='none')
    ax.axhline(threshold, color='#d62728', linestyle='--', linewidth=1.5)
    ax.axvline(0, color='gray', linewidth=1.2)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('log2FC', fontsize=11)
    ax.set_ylabel('-log10 p-value', fontsize=11)
    ax.grid(True, alpha=0.4)

# panel-specific visible ranges approximated from screenshot
axes[0].set_xlim(-3.8, 4.1); axes[0].set_ylim(-0.2, 4.1)
axes[1].set_xlim(-6.2, 6.2); axes[1].set_ylim(-0.5, 10.5)
axes[2].set_xlim(-1.9, 4.3); axes[2].set_ylim(-0.2, 3.7)
axes[3].set_xlim(-3.2, 6.4); axes[3].set_ylim(-0.5, 10.5)
axes[4].set_xlim(-3.8, 1.8); axes[4].set_ylim(-0.1, 2.65)
axes[5].set_xlim(-1.55, 2.35); axes[5].set_ylim(-0.1, 3.05)

plt.tight_layout()
plt.savefig('organelle_volcano_panels_embedded.png', bbox_inches='tight')
plt.savefig('organelle_volcano_panels_embedded.pdf', bbox_inches='tight')
plt.show()

for i, (title, df) in enumerate(panels.items(), start=1):
    safe = title.split(' (')[0].lower().replace(' ', '_')
    df.to_csv(f'{safe}_panel_points.csv', index=False)

print('Saved files: organelle_volcano_panels_embedded.png, .pdf, and per-panel CSVs')
