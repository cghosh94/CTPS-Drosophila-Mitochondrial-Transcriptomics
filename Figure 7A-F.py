# ======================================================================
# Gene-Drug Interaction Figures (TIFF) - Upregulated & Downregulated
# CORRECTED: Removed "Antioxidant" from MDH2 interactions
# Outputs: Heatmap, Dotplot, Network, Barplot for each group
# ======================================================================

!pip install networkx -q
!pip install matplotlib seaborn pandas networkx -q

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ----------------------------------------------------------------
# 1. ORIGINAL GENE DATA (Log2FC and p-values from your list)
# ----------------------------------------------------------------

# Upregulated genes
up_genes_original = [
    "CRAT", "MetRS-m", "ScpX", "IscU", "GatB", "Tom7", "Trap1",
    "ND-13A", "Fis1", "ND-B14", "Hsc20", "COX7AL", "AlkB", "GluRS-m",
    "Lon", "Fdx1", "Bmcp", "ArgRS-m", "mtRNApol", "Atg3", "mtSSB",
    "D2hgdh", "ND-13B", "Suv3", "Gdh", "Tom20", "Ttc19", "Alr", "Letm1"
]
up_log2fc = [
    0.718, 0.416, 0.448, 0.832, 0.610, 0.190, 0.360,
    0.139, 0.260, 0.213, 0.536, 1.766, 0.667, 0.806,
    0.799, 0.576, 0.341, 0.204, 0.254, 0.190, 0.082,
    0.517, 0.266, 0.224, 0.518, 0.129, 0.646, 0.718, 0.264
]
up_pvals = [
    0.00029874, 0.0022835, 0.0027442, 0.0035896, 0.0039404, 0.0067562, 0.010176,
    0.010615, 0.011611, 0.012505, 0.013162, 0.013236, 0.014035, 0.014755,
    0.015515, 0.015574, 0.018589, 0.020736, 0.022538, 0.027217, 0.02928,
    0.032235, 0.035119, 0.035136, 0.03599, 0.04089, 0.046433, 0.046479, 0.049017
]
up_gene_dict = dict(zip(up_genes_original, up_log2fc))

# Downregulated genes
down_genes_original = [
    "sesB", "QIL1", "Mpcp2", "ATPsynG", "iPLA2-VIA", "aralar1",
    "Mdh2", "Hsp60A", "Etf-QO", "mEFTu1", "ND-ACP", "blw", "ND-30"
]
down_log2fc = [
    -0.361, -0.211, -0.375, -0.157, -0.615, -0.635,
    -0.299, -0.353, -0.416, -0.239, -0.201, -0.171, -0.222
]
down_pvals = [
    0.002528, 0.0064847, 0.0077837, 0.0080681, 0.010421, 0.011768,
    0.015328, 0.018416, 0.019111, 0.026075, 0.032309, 0.039089, 0.042545
]
down_gene_dict = dict(zip(down_genes_original, down_log2fc))

# ----------------------------------------------------------------
# 2. INTERACTION DICTIONARY (Gene -> list of (Drug, Approval Status))
# CORRECTED: Removed "ANTIOXIDANT" from Mdh2
# ----------------------------------------------------------------

interactions = {
    # Upregulated
    "CRAT":       [("Levocarnitine", "Approved"), ("Levocarnitine Propionate", "Approved")],
    "MetRS-m":    [("MetRS inhibitor", "Not Approved")],
    "IscU":       [("Erlotinib", "Approved")],
    "Trap1":      [("TRAP1-IN-1", "Not Approved")],
    "Fis1":       [("P110", "Not Approved")],
    "AlkB":       [("Alectinib", "Approved")],
    "GluRS-m":    [("GluRS inhibitor", "Not Approved")],
    "Lon":        [("Lonafarnib", "Approved")],
    "Fdx1":       [
        ("TAS05567", "Approved"), ("LANRAPLENIB", "Approved"),
        ("SYK INHIBITOR", "Approved"), ("ER-27319", "Approved"),
        ("GUSACITINIB", "Approved"), ("GSK143", "Approved"),
        ("SYK INHIBITOR II", "Approved"), ("SOVLEPLENIB", "Approved"),
        ("TAMATINIB", "Approved"), ("ENTOSPLETINIB", "Approved"),
        ("COMPOUND 7", "Approved"), ("COMPOUND 23", "Approved"),
        ("MIVAVOTINIB", "Approved"), ("LAZERTINIB", "Approved"),
        ("P505-15", "Approved"), ("GSK2646264", "Approved"),
        ("COMPOUND 12D", "Approved")
    ],
    "Bmcp":       [("Trichostatin A", "Not Approved")],
    "ArgRS-m":    [("Borrelidin", "Not Approved")],
    "mtRNApol":   [("Azacytidine", "Approved")],
    "Atg3":       [("Chloroquine", "Approved")],
    "D2hgdh":     [("Enasidenib", "Approved")],
    "Gdh":        [("Purpurin", "Not Approved")],
    "Alr":        [("Ponalrestat", "Not Approved")],
    
    # Downregulated (CORRECTED: ANTIOXIDANT removed)
    "ATPsynG":    [("Crizotinib", "Approved")],
    "iPLA2-VIA":  [("FKGK18", "Not Approved")],
    "Mdh2":       [
        ("EXPECTORANT", "Not Approved"),
        ("THERAPEUTIC GONADOTROPIN", "Approved"),
        ("NIFLUMIC ACID", "Approved"),
        ("ETOPOSIDE", "Approved"),
        ("CISPLATIN", "Approved")
        # "ANTIOXIDANT" removed (not a suitable drug)
    ],
    "Hsp60A":     [("Mizoribine", "Approved")],
    "Etf-QO":     [("Riboflavin", "Approved")],
    "blw":        [("Navitoclax", "Not Approved")],
}

# ----------------------------------------------------------------
# 3. BUILD DATAFRAMES
# ----------------------------------------------------------------

def build_interaction_df(gene_list, log2fc_dict, interactions):
    rows = []
    for gene in gene_list:
        if gene in interactions and interactions[gene]:
            for drug, status in interactions[gene]:
                rows.append({
                    'Gene': gene,
                    'Drug': drug,
                    'Status': status,
                    'Log2FC': log2fc_dict.get(gene, 0)
                })
        else:
            # Gene with no known drug – keep for completeness but filter later
            rows.append({
                'Gene': gene,
                'Drug': 'None',
                'Status': 'Unknown',
                'Log2FC': log2fc_dict.get(gene, 0)
            })
    return pd.DataFrame(rows)

df_up_all = build_interaction_df(up_genes_original, up_gene_dict, interactions)
df_down_all = build_interaction_df(down_genes_original, down_gene_dict, interactions)

# Keep only genes with at least one real drug
df_up = df_up_all[df_up_all['Drug'] != 'None'].copy()
df_down = df_down_all[df_down_all['Drug'] != 'None'].copy()

print("✅ Upregulated interactions:", len(df_up))
print("✅ Downregulated interactions:", len(df_down))
print("\nUpregulated genes with drugs:", df_up['Gene'].unique())
print("Downregulated genes with drugs:", df_down['Gene'].unique())

# ----------------------------------------------------------------
# 4. PLOTTING FUNCTIONS
# ----------------------------------------------------------------

def plot_heatmap(df, title, filename, is_up=True):
    """Binary heatmap: genes x drugs, with color bars for status."""
    if df.empty:
        print(f"⚠️ No data for {title}")
        return
    
    # Pivot to matrix
    matrix = df.pivot_table(index='Gene', columns='Drug', values='Log2FC', aggfunc='count', fill_value=0)
    matrix = (matrix > 0).astype(int)
    
    # Row colors (regulation)
    row_color = '#ff7f0e' if is_up else '#1f77b4'
    row_colors = pd.Series([row_color]*len(matrix), index=matrix.index, name="Regulation")
    
    # Column colors by approval status
    drug_status = df.drop_duplicates('Drug').set_index('Drug')['Status']
    col_colors = pd.Series(
        ['#2ca02c' if drug_status.get(d, 'Unknown') == 'Approved' else '#d62728' for d in matrix.columns],
        index=matrix.columns,
        name="Approval"
    )
    
    figsize = (max(10, len(matrix.columns)*0.4), max(6, len(matrix)*0.4))
    g = sns.clustermap(
        matrix,
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        col_colors=col_colors,
        cmap='Blues',
        linewidths=0.5,
        linecolor='lightgrey',
        cbar_kws={'label': 'Interaction (1=Yes, 0=No)', 'shrink': 0.6},
        figsize=figsize,
        annot=True,
        fmt='d',
        annot_kws={'size': 8}
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=9)
    g.ax_heatmap.set_title(title, pad=20, fontsize=14)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, format='tiff', bbox_inches='tight')
    plt.show()
    print(f"✅ Heatmap saved: {filename}")

def plot_dotplot(df, title, filename):
    """Dot plot: Log2FC on x, gene on y, each drug as point colored by status."""
    if df.empty:
        print(f"⚠️ No data for {title}")
        return
    
    genes = df['Gene'].unique()
    y_pos = {gene: i for i, gene in enumerate(genes)}
    df['y'] = df['Gene'].map(y_pos)
    # Jitter to separate drugs per gene
    df['y_jitter'] = df['y'] + (df.groupby('Gene').cumcount() - (df.groupby('Gene')['Gene'].transform('size') - 1)/2) * 0.12
    
    fig, ax = plt.subplots(figsize=(12, max(6, len(genes)*0.5)))
    
    color_map = {'Approved': '#2ca02c', 'Not Approved': '#d62728', 'Unknown': '#7f7f7f'}
    colors = df['Status'].map(color_map)
    
    scatter = ax.scatter(df['Log2FC'], df['y_jitter'], c=colors, s=100, alpha=0.85, edgecolors='black', linewidth=0.5)
    
    # Add drug labels
    for _, row in df.iterrows():
        ax.text(row['Log2FC'] + 0.02, row['y_jitter'], row['Drug'], fontsize=7, va='center', alpha=0.8)
    
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(list(y_pos.keys()), fontsize=10)
    ax.set_xlabel('Log2 Fold Change (CTPs_Ri vs Control)', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    ax.set_title(title, fontsize=14)
    
    handles = [Patch(color=c, label=s) for s, c in color_map.items() if s in df['Status'].unique()]
    ax.legend(handles=handles, title='Drug Status', loc='lower right')
    ax.grid(axis='x', linestyle=':', alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, format='tiff', bbox_inches='tight')
    plt.show()
    print(f"✅ Dotplot saved: {filename}")

def plot_network(df, title, filename, is_up=True):
    """Circular network: genes and drugs as nodes."""
    if df.empty:
        print(f"⚠️ No data for {title}")
        return
    
    G = nx.Graph()
    genes = df['Gene'].unique()
    drugs = df['Drug'].unique()
    drug_status = df.drop_duplicates('Drug').set_index('Drug')['Status'].to_dict()
    
    for g in genes:
        G.add_node(g, type='gene')
    for d in drugs:
        G.add_node(d, type='drug', status=drug_status.get(d, 'Unknown'))
    
    for _, row in df.iterrows():
        G.add_edge(row['Gene'], row['Drug'])
    
    pos = nx.circular_layout(G)
    
    gene_nodes = [n for n, attr in G.nodes(data=True) if attr['type']=='gene']
    drug_nodes = [n for n, attr in G.nodes(data=True) if attr['type']=='drug']
    
    drug_colors = []
    for d in drug_nodes:
        status = G.nodes[d]['status']
        if status == 'Approved':
            drug_colors.append('#2ca02c')
        elif status == 'Not Approved':
            drug_colors.append('#d62728')
        else:
            drug_colors.append('#7f7f7f')
    
    fig, ax = plt.subplots(figsize=(12, 12))
    
    gene_color = '#ff7f0e' if is_up else '#1f77b4'
    nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_color=gene_color, node_size=800, alpha=0.9)
    nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes, node_color=drug_colors, node_size=600, alpha=0.8)
    nx.draw_networkx_edges(G, pos, width=1.5, alpha=0.5, edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
    
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=gene_color, markersize=10, label='Genes'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ca02c', markersize=10, label='Approved Drugs'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728', markersize=10, label='Not Approved Drugs'),
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=10)
    ax.set_title(title, fontsize=16)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, format='tiff', bbox_inches='tight')
    plt.show()
    print(f"✅ Network saved: {filename}")

def plot_barplot(df, title, filename):
    """Stacked bar: number of drugs per gene, split by approval status."""
    if df.empty:
        print(f"⚠️ No data for {title}")
        return
    
    counts = df.groupby(['Gene', 'Status']).size().unstack(fill_value=0)
    counts = counts.reindex(df['Gene'].unique(), fill_value=0)
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(counts)*0.4)))
    counts.plot(kind='barh', stacked=True, ax=ax, color={'Approved':'#2ca02c', 'Not Approved':'#d62728', 'Unknown':'#7f7f7f'})
    ax.set_xlabel('Number of interacting drugs', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(title='Drug Status')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, format='tiff', bbox_inches='tight')
    plt.show()
    print(f"✅ Barplot saved: {filename}")

# ----------------------------------------------------------------
# 5. GENERATE FIGURES FOR UPREGULATED
# ----------------------------------------------------------------
print("\n" + "="*60)
print("GENERATING UPREGULATED FIGURES")
print("="*60)
if not df_up.empty:
    plot_heatmap(df_up, "Upregulated Genes – Drug Interaction Matrix", "Up_Heatmap.tiff", is_up=True)
    plot_dotplot(df_up, "Upregulated Genes – Log2FC vs Drugs", "Up_Dotplot.tiff")
    plot_network(df_up, "Upregulated Genes – Drug Network", "Up_Network.tiff", is_up=True)
    plot_barplot(df_up, "Upregulated Genes – Drug Count per Gene", "Up_Barplot.tiff")
else:
    print("⚠️ No upregulated interactions to plot.")

# ----------------------------------------------------------------
# 6. GENERATE FIGURES FOR DOWNREGULATED
# ----------------------------------------------------------------
print("\n" + "="*60)
print("GENERATING DOWNREGULATED FIGURES")
print("="*60)
if not df_down.empty:
    plot_heatmap(df_down, "Downregulated Genes – Drug Interaction Matrix", "Down_Heatmap.tiff", is_up=False)
    plot_dotplot(df_down, "Downregulated Genes – Log2FC vs Drugs", "Down_Dotplot.tiff")
    plot_network(df_down, "Downregulated Genes – Drug Network", "Down_Network.tiff", is_up=False)
    plot_barplot(df_down, "Downregulated Genes – Drug Count per Gene", "Down_Barplot.tiff")
else:
    print("⚠️ No downregulated interactions to plot.")

print("\n" + "="*60)
print("🎉 ALL TIFF FIGURES GENERATED SUCCESSFULLY!")
print("="*60)
print("\n📁 Files generated:")
print("   Up_Heatmap.tiff, Up_Dotplot.tiff, Up_Network.tiff, Up_Barplot.tiff")
print("   Down_Heatmap.tiff, Down_Dotplot.tiff, Down_Network.tiff, Down_Barplot.tiff")
print("\n📥 Download from Colab File pane (left sidebar) → Right-click → Download")
