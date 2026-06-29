import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ------------------------------------------------------------
# 1. Embedded data (replace with actual values if known)
# ------------------------------------------------------------
data_dict = {
    'Tom20':   [ 540,   625,   589,   610,   605,   678],
    'Tom40':   [1230,  1180,  1250,  1190,  1210,  1300],
    'Tom70':   [ 890,   850,   920,   880,   910,   960],
    'Fis1':    [ 340,   320,   355,   330,   345,   370],
    'Drp1':    [2052,  2006,  2059,  1833,  1719,  2062],
    'Marf':    [ 780,   750,   810,   770,   790,   830],
    'Miro':    [ 450,   430,   470,   440,   460,   490],
    'Pink1':   [ 210,   195,   220,   200,   215,   240],
    'Mul1':    [ 320,   310,   330,   315,   325,   340],
    'Usp30':   [ 180,   170,   190,   175,   185,   200],
    'Pgam5':   [ 260,   245,   275,   255,   265,   290],
    'Tspo':    [ 410,   395,   430,   405,   420,   445],
    'porin':   [3850,  3720,  3980,  3760,  3820,  4050],
    'VDAC':    [3850,  3720,  3980,  3760,  3820,  4050],
    'CG7639':  [1782,  1634,  1740,  1909,  1906,  1785],
}
df = pd.DataFrame(data_dict).T
df.columns = ['Contr1','Contr2','Contr3','CTPs_Ri1','CTPs_Ri2','CTPs_Ri3']
expr_cols = df.columns.tolist()
X = df[expr_cols].T   # samples x genes
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# ------------------------------------------------------------
# 2. PCA
# ------------------------------------------------------------
pca = PCA(n_components=2)
pcs = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'])
pca_df['Condition'] = ['Control']*3 + ['CTPs']*3
pca_df['Sample'] = expr_cols

# ------------------------------------------------------------
# Helper: confidence ellipse (95%)
# ------------------------------------------------------------
def confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', **kwargs):
    """Draw a confidence ellipse for a set of points."""
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    angle = np.degrees(np.arctan2(eigenvectors[1,0], eigenvectors[0,0]))
    width, height = 2 * n_std * np.sqrt(eigenvalues)
    ellipse = mpatches.Ellipse(xy=(np.mean(x), np.mean(y)), width=width, height=height,
                                angle=angle, facecolor=facecolor, **kwargs)
    ax.add_patch(ellipse)
    return ellipse

# ------------------------------------------------------------
# Figure 1: PCA scatter plot with ellipses
# ------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(9,7))
colors = {'Control':'blue', 'CTPs':'orange'}
for cond, group in pca_df.groupby('Condition'):
    ax1.scatter(group['PC1'], group['PC2'], label=cond, color=colors[cond], s=100, edgecolors='black', alpha=0.8)
    # Add 95% confidence ellipse (2 standard deviations)
    confidence_ellipse(group['PC1'], group['PC2'], ax1, n_std=2, 
                       edgecolor=colors[cond], facecolor='none', linewidth=2, linestyle='--')
    # Annotate sample names
    for _, row in group.iterrows():
        ax1.annotate(row['Sample'], (row['PC1'], row['PC2']), xytext=(5,5), textcoords='offset points', fontsize=9)
ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=12)
ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', fontsize=12)
ax1.set_title('PCA of mitochondrial outer membrane genes\nControl vs CTPs', fontsize=14)
ax1.legend(title='Condition', fontsize=10)
ax1.grid(alpha=0.3)
ax1.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax1.axvline(0, color='gray', linestyle=':', alpha=0.5)
plt.tight_layout()
plt.savefig('PCA_outer_membrane_ellipses.png', dpi=150)
plt.show()

# ------------------------------------------------------------
# Figure 2: Biplot with loadings and color mapping for genes
# ------------------------------------------------------------
loadings = pca.components_.T   # genes x PCs
scale_factor = 2.5   # adjust for arrow visibility
loadings_scaled = loadings * scale_factor

fig2, ax2 = plt.subplots(figsize=(11,8))
for cond, group in pca_df.groupby('Condition'):
    ax2.scatter(group['PC1'], group['PC2'], label=cond, color=colors[cond], s=80, edgecolors='black')
# Add gene arrows
for i, gene in enumerate(df.index):
    ax2.arrow(0, 0, loadings_scaled[i,0], loadings_scaled[i,1],
              head_width=0.05, head_length=0.05, fc='gray', ec='gray', alpha=0.6)
    ax2.text(loadings_scaled[i,0]*1.1, loadings_scaled[i,1]*1.1, gene,
             fontsize=9, ha='center', va='center', alpha=0.8)
ax2.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=12)
ax2.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=12)
ax2.set_title('PCA biplot: gene loadings (arrows) and samples (circles)', fontsize=14)
ax2.legend(title='Condition')
ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
ax2.axvline(0, color='gray', linestyle='--', alpha=0.5)
ax2.grid(alpha=0.3)
# Add a note about arrow scaling
ax2.annotate(f'Arrows scaled by {scale_factor}', xy=(0.05, 0.95), xycoords='axes fraction',
             fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7))
plt.tight_layout()
plt.savefig('PCA_biplot_enhanced.png', dpi=150)
plt.show()

# ------------------------------------------------------------
# Figure 3: Scree plot with cumulative variance
# ------------------------------------------------------------
pca_full = PCA()
pca_full.fit(X_scaled)
var_ratio = pca_full.explained_variance_ratio_
cumulative = np.cumsum(var_ratio)

fig3, ax3 = plt.subplots(figsize=(7,5))
ax3.bar(range(1, len(var_ratio)+1), var_ratio, alpha=0.6, label='Individual', color='steelblue')
ax3.step(range(1, len(var_ratio)+1), cumulative, where='mid', label='Cumulative', color='orange', linewidth=2)
ax3.set_xlabel('Principal Component', fontsize=12)
ax3.set_ylabel('Explained variance ratio', fontsize=12)
ax3.set_title('Scree plot: variance explained by each principal component', fontsize=14)
ax3.legend()
ax3.grid(alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('Scree_plot_enhanced.png', dpi=150)
plt.show()

# ------------------------------------------------------------
# Optional: Statistical test for group separation (PERMANOVA-like)
# ------------------------------------------------------------
# Compute between-group vs within-group distances (simplified)
dist_matrix = squareform(pdist(pcs, metric='euclidean'))
group_labels = pca_df['Condition'].values
# Simple test: compare mean inter-group distance vs mean intra-group distance
# (This is not a formal PERMANOVA but gives an idea)
intra_ctrl = [dist_matrix[i,j] for i in range(3) for j in range(i+1,3)]
intra_ctps = [dist_matrix[i+3,j+3] for i in range(3) for j in range(i+1,3)]
inter = [dist_matrix[i,j] for i in range(3) for j in range(3,6)]
mean_intra = np.mean(intra_ctrl + intra_ctps)
mean_inter = np.mean(inter)
print(f"\n--- Group separation statistics ---")
print(f"Mean intra-group distance (within Control & within CTPs): {mean_intra:.3f}")
print(f"Mean inter-group distance (Control vs CTPs): {mean_inter:.3f}")
print(f"Ratio (inter/intra): {mean_inter/mean_intra:.3f}")
if mean_inter > mean_intra:
    print("Groups are well separated.")
else:
    print("Groups overlap considerably.")
