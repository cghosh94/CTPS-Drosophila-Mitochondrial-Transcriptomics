#!/usr/bin/env python3
# ============================================================
# VOLCANO PLOT: COLOR-CODED BY REGULATION (UP/DOWN)
# - Red: up-regulated significant (log2FC > +1, p < 0.05)
# - Blue: down-regulated significant (log2FC < -1, p < 0.05)
# - Grey: non-significant or |log2FC| <= 1
# Labels all significant genes (both directions)
# No external file – embedded realistic gene symbols
# ============================================================

import subprocess
import sys
import importlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ------------------ Auto-install adjustText if missing ------------------
def install_and_import(package):
    try:
        importlib.import_module(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    finally:
        return importlib.import_module(package)

adjustText = install_and_import("adjustText")
from adjustText import adjust_text
# -----------------------------------------------------------------------

# ================== USER PARAMETERS ==================
OUTPUT_PREFIX = "volcano_up_down_colored"
P_THRESHOLD = 0.05
FC_THRESHOLD = 1.0          # |log2FC| > 1 for significance
RANDOM_SEED = 42
N_GENES = 200
# =====================================================

np.random.seed(RANDOM_SEED)

# ------------------ Realistic gene symbol list (CG codes + common names) ------------------
real_gene_symbols = [
    "CG31683", "CG7024", "CG6012", "CG14806", "CG5044", "CG5776", "CG4908", "CG4942",
    "CG4586", "CG14812", "CG5214", "CG5028", "CG42376", "CG33182", "CG1542", "CG8878",
    "CG10684", "CG17278", "CG18179", "CG11073", "CG13096", "CG18311", "CG18617", "CG12262",
    "ND-51L2", "sad", "sro", "COX7AL", "fs(1)Ya", "Qtzl", "Cyp12d1-d", "ACC", "CPT2",
    "Acsl", "bgm", "Cyt-c", "IMPDH2", "HIBCH", "ATPsyn-β", "COX5A", "Drp1", "Acon",
    "COX4", "SPATA5", "BCS1", "COX18", "Acyl-CoA oxidase", "Ank2", "IMPDH", "DLST",
    "IDHy", "CTPsyn", "ND-42", "ND-75", "ND-49", "ND-20", "ND-18", "ND-24", "ND-15",
    "ND-ML1", "ND-ML2", "ND-51L1", "SDHA", "SDHB", "SDHC", "SDHD", "UQCR-C1",
    "UQCR-Q", "UQCR-14", "COX1", "COX2", "COX3", "COX4L", "COX5B", "COX6A", "COX6B",
    "COX7A", "COX7B", "COX8", "ATP5A", "ATP5B", "ATP5C", "ATP5D", "ATP5E", "ATP5F",
    "ATP5G", "ATP5H", "ATP5J", "ATP5L", "ANT", "TCA", "IDH3A", "IDH3B", "IDH3G", "OGDH",
    "DLD", "PDHA1", "PDHB", "PDK1", "PDK2", "PDK3", "PDK4", "GOT1", "GOT2", "MDH1", "MDH2"
]

# Ensure enough unique names
if len(real_gene_symbols) < N_GENES:
    extended = real_gene_symbols.copy()
    while len(extended) < N_GENES:
        for sym in real_gene_symbols:
            extended.append(f"{sym}_dup")
            if len(extended) >= N_GENES:
                break
    gene_symbols = extended[:N_GENES]
else:
    gene_symbols = real_gene_symbols[:N_GENES]

np.random.shuffle(gene_symbols)

# ------------------ Generate synthetic expression data ------------------
log2FC = np.concatenate([
    np.random.normal(0, 0.5, int(N_GENES * 0.7)),          # background
    np.random.normal(2.5, 0.8, int(N_GENES * 0.15)),       # up-regulated
    np.random.normal(-2.3, 0.7, int(N_GENES * 0.15))       # down-regulated
])[:N_GENES]

p_value = np.zeros(N_GENES)
for i in range(N_GENES):
    if abs(log2FC[i]) > 1.5:
        p_value[i] = np.random.exponential(0.01)
    elif abs(log2FC[i]) > 0.8:
        p_value[i] = np.random.uniform(0.01, 0.2)
    else:
        p_value[i] = np.random.uniform(0.1, 0.9)

# Extra significant hits
p_value[np.abs(log2FC) > 2] = np.random.uniform(1e-6, 0.001, sum(np.abs(log2FC) > 2))

df = pd.DataFrame({
    'gene': gene_symbols,
    'log2FC': log2FC,
    'p_value': p_value
})
df['neglog10_pvalue'] = -np.log10(df['p_value'])

# ------------------ Classify significance and direction ------------------
sig_up = df[(df['p_value'] < P_THRESHOLD) & (df['log2FC'] > FC_THRESHOLD)]
sig_down = df[(df['p_value'] < P_THRESHOLD) & (df['log2FC'] < -FC_THRESHOLD)]
non_sig = df[~((df['p_value'] < P_THRESHOLD) & (np.abs(df['log2FC']) > FC_THRESHOLD))]

print(f"Total genes: {len(df)}")
print(f"Significant up-regulated (p<{P_THRESHOLD}, log2FC>{FC_THRESHOLD}): {len(sig_up)}")
print(f"Significant down-regulated (p<{P_THRESHOLD}, log2FC<{-FC_THRESHOLD}): {len(sig_down)}")

# ------------------ Create volcano plot ------------------
plt.figure(figsize=(12, 8))

# Non-significant (grey)
plt.scatter(non_sig['log2FC'], non_sig['neglog10_pvalue'],
            c='lightgrey', s=30, alpha=0.6, label='Not significant')

# Significant up (red)
if len(sig_up) > 0:
    plt.scatter(sig_up['log2FC'], sig_up['neglog10_pvalue'],
                c='red', s=50, alpha=0.8, label='Up-regulated (sig)')

# Significant down (blue)
if len(sig_down) > 0:
    plt.scatter(sig_down['log2FC'], sig_down['neglog10_pvalue'],
                c='blue', s=50, alpha=0.8, label='Down-regulated (sig)')

# Threshold lines
plt.axhline(-np.log10(P_THRESHOLD), color='black', linestyle='--',
            linewidth=1.5, label=f'p = {P_THRESHOLD}')
plt.axvline(0, color='gray', linestyle='--', linewidth=1)
plt.axvline(-FC_THRESHOLD, color='gray', linestyle=':', linewidth=1)
plt.axvline(FC_THRESHOLD, color='gray', linestyle=':', linewidth=1)

# ------------------ Label all significant genes (both up and down) ------------------
texts = []
sig_all = pd.concat([sig_up, sig_down])
for _, row in sig_all.iterrows():
    txt = plt.text(row['log2FC'], row['neglog10_pvalue'],
                   row['gene'], fontsize=8, ha='center', va='center')
    texts.append(txt)

if len(texts) > 0:
    adjust_text(texts,
                arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
                expand_points=(1.5, 1.5),
                expand_text=(1.2, 1.2))

# ------------------ Format and save ------------------
plt.xlabel('log2 Fold Change', fontsize=12)
plt.ylabel('-log10(p-value)', fontsize=12)
plt.title(f'Volcano plot: up-regulated (red) vs down-regulated (blue)\n(p < {P_THRESHOLD}, |log2FC| > {FC_THRESHOLD})', fontsize=12)
plt.legend(loc='upper right')
plt.tight_layout()

plt.savefig(f"{OUTPUT_PREFIX}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTPUT_PREFIX}.pdf", bbox_inches='tight')
plt.show()

# Save significant genes (with direction)
sig_all.to_csv(f"{OUTPUT_PREFIX}_significant_genes.csv", index=False)
print(f"\nSaved files:\n - {OUTPUT_PREFIX}.png\n - {OUTPUT_PREFIX}.pdf\n - {OUTPUT_PREFIX}_significant_genes.csv")
