import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.lines import Line2D

# 3 red samples + 3 blue samples
ctps = np.array([1400, 1450, 2200, 3400, 3650, 3780])
groups = np.array(['red', 'red', 'red', 'blue', 'blue', 'blue'])

panel_data = {
    'ACC': np.array([60500, 59000, 101000, 115500, 157000, 140000]),
    'CPT2': np.array([1505, 1510, 1755, 1860, 1920, 1850]),
    'Acsl': np.array([19400, 17500, 27300, 30000, 32100, 31600]),
    'bgm': np.array([6450, 5750, 11750, 12350, 16050, 16200]),
    'Cyt-c': np.array([23000, 23100, 23600, 26500, 30200, 28900]),
    'CRAT': np.array([1640, 1530, 1510, 1000, 960, 890]),
    'Catalase': np.array([4030, 3230, 2710, 2550, 1960, 2210]),
    'ECI (CG4592)': np.array([1225, 1090, 1060, 750, 950, 800])
}

shown_stats = {
    'ACC': ('0.964', '1.95e-03'),
    'CPT2': ('0.957', '2.70e-03'),
    'Acsl': ('0.953', '3.22e-03'),
    'bgm': ('0.953', '3.31e-03'),
    'Cyt-c': ('0.945', '4.43e-03'),
    'CRAT': ('-0.979', '6.79e-04'),
    'Catalase': ('-0.902', '1.38e-02'),
    'ECI (CG4592)': ('-0.869', '2.48e-02')
}

ylims_ticks = {
    'ACC': ((55000, 162000), [60000, 80000, 100000, 120000, 140000, 160000]),
    'CPT2': ((1485, 1940), [1500, 1600, 1700, 1800, 1900]),
    'Acsl': ((16800, 33000), [17500, 20000, 22500, 25000, 27500, 30000, 32500]),
    'bgm': ((5200, 16600), [6000, 8000, 10000, 12000, 14000, 16000]),
    'Cyt-c': ((22200, 30600), [24000, 26000, 28000, 30000]),
    'CRAT': ((850, 1670), [1000, 1200, 1400, 1600]),
    'Catalase': ((1900, 4100), [2000, 2500, 3000, 3500, 4000]),
    'ECI (CG4592)': ((720, 1250), [800, 900, 1000, 1100, 1200])
}

x_ticks = [1500, 2000, 2500, 3000, 3500]

fig, axes = plt.subplots(2, 4, figsize=(15.36, 7.68))
fig.suptitle('Correlation of CTPS with Key Mitochondrial Genes', fontsize=20, y=0.98)
axes = axes.flatten()

for ax, (gene, y) in zip(axes, panel_data.items()):
    red_mask = groups == 'red'
    blue_mask = groups == 'blue'

    ax.scatter(
        ctps[red_mask], y[red_mask],
        color='#ff3b3b', edgecolor='black', s=80, alpha=0.8, label='Red group'
    )
    ax.scatter(
        ctps[blue_mask], y[blue_mask],
        color='#3b44ff', edgecolor='black', s=80, alpha=0.8, label='Blue group'
    )

    slope, intercept, r, p, _ = linregress(ctps, y)
    xline = np.linspace(1400, 3800, 200)
    ax.plot(xline, slope * xline + intercept, linestyle='--', color='gray', linewidth=1.6)

    r_txt, p_txt = shown_stats[gene]
    ax.set_title(f'{gene}\nr = {r_txt}, p = {p_txt}', fontsize=17)
    ax.set_xlabel('CTPS expression (counts)', fontsize=14)

    if gene == 'ACC':
        ylabel = 'ACC expression'
    elif gene == 'Cyt-c':
        ylabel = 'Cyt-c expression'
    elif gene == 'Catalase':
        ylabel = 'Catalase expression'
    elif gene == 'bgm':
        ylabel = 'bgm expression'
    elif gene == 'CPT2':
        ylabel = 'CPT2 expression'
    elif gene == 'Acsl':
        ylabel = 'Acsl expression'
    elif gene == 'CRAT':
        ylabel = 'CRAT expression'
    elif gene == 'ECI (CG4592)':
        ylabel = 'ECI (CG4592) expression'
    else:
        ylabel = gene

    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlim(1300, 3900)
    ax.set_xticks(x_ticks)

    ylim, yticks = ylims_ticks[gene]
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.grid(True, alpha=0.7)

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Red group',
           markerfacecolor='#ff3b3b', markeredgecolor='black', markersize=9),
    Line2D([0], [0], marker='o', color='w', label='Blue group',
           markerfacecolor='#3b44ff', markeredgecolor='black', markersize=9),
    Line2D([0], [0], linestyle='--', color='gray', label='Linear fit')
]

fig.legend(
    handles=legend_elements,
    loc='upper center',
    ncol=3,
    frameon=False,
    bbox_to_anchor=(0.5, 1.02),
    fontsize=12
)

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.show()
