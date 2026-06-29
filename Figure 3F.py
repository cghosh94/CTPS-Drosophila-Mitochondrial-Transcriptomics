import numpy as np
import matplotlib.pyplot as plt

# Embedded numerical data reconstructed from the radar plot
labels = [
    'ACC', 'CPT2', 'Acsl', 'bgm', 'Cyt-c', 'IMPDH2 (CG14806)',
    'HIBCH (CG5044)', 'ATPsyn-β', 'CRAT', 'Catalase', 'ECI (CG4592)',
    'FUNDC1 (CG5676)', 'Atg8a', 'cora (protein 4.1)', 'MTO1 (CG4610)', 'ClpP (CG5045)'
]

values = [
    0.96, 0.95, 0.94, 0.95, 0.97, 0.82,
    0.78, 0.72, -0.98, -0.95, -0.92,
    -0.88, -0.82, -0.76, -0.62, -0.45
]

N = len(labels)
angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
angles += angles[:1]
values_closed = values + values[:1]

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
ax.plot(angles, values_closed, linewidth=2.2, marker='o')
ax.fill(angles, values_closed, alpha=0.25)

ax.set_xticks(angles[:-1])
ax.set_xticklabels(labels, fontsize=10)
ax.set_yticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticklabels(['-1', '-0.5', '0', '0.5', '1'])
ax.set_ylim(-1, 1)
ax.set_title('Radar Plot: Top Positive and Negative Correlations with CTPS', fontsize=16)

plt.tight_layout()
plt.show()
