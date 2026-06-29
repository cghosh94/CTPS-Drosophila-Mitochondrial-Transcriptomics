CTPS-Mitochondrial-Transcriptome
Transcriptomic Re-analysis Reveals Coordinated Mitochondrial Remodeling Following CTPS Depletion in Drosophila melanogaster

This repository contains the custom Python scripts and associated resources used for the computational analyses presented in our manuscript:

"Transcriptomic Re-analysis Reveals Coordinated Remodeling of the Mitochondrial Transcriptome in CTPS-Depleted Drosophila Fat Body."

The study systematically re-analyzes publicly available RNA-sequencing data to investigate the impact of CTP synthase (CTPS) depletion on mitochondrial transcriptional programs. Comprehensive bioinformatic analyses revealed coordinated alterations in mitochondrial metabolism, fatty acid oxidation, oxidative phosphorylation, mitochondrial translation, outer mitochondrial membrane organization, and pharmacologically targetable mitochondrial networks.

Dataset

Source: NCBI Gene Expression Omnibus (GEO)

Accession: GSE221707

Organism: Drosophila melanogaster

Tissue: Fat body

Experimental groups:

Control
CTPS-RNAi

Input data used in this study:

Raw gene count matrix (counts annotation file)
Sample annotation file

Note: Differential expression and all downstream analyses were performed using the raw count matrix. The normalized FPKM file provided with the dataset was not used.

Analyses Included

The repository contains Python scripts for:

Differential gene expression analysis
Organelle-specific transcriptomic profiling
Principal component analysis (PCA)
Pearson correlation analysis
Linear regression analysis
Heatmap generation
Volcano plot visualization
KEGG pathway enrichment analysis
Gene Ontology (GO) enrichment analysis
Reactome pathway analysis
STRING protein–protein interaction analysis
Drug–gene interaction analysis
Publication-quality figure generation
Software Requirements

The analyses were performed using Python in Google Colab.

Major Python packages include:

NumPy
Pandas
SciPy
Matplotlib
Seaborn
Scikit-learn
Plotly
NetworkX
Statsmodels

Additional package requirements are listed in requirements.txt.

Reproducibility

All scripts are provided to facilitate reproducibility of the computational analyses and publication figures. Unless otherwise indicated within individual scripts, analyses were performed using the raw count matrix from GSE221707 together with the corresponding sample annotation file.

Citation

If you use this repository, please cite:

Ghosh C., et al. Transcriptomic Re-analysis Reveals Coordinated Remodeling of the Mitochondrial Transcriptome in CTPS-Depleted Drosophila Fat Body. (Manuscript under review.)

Contact

Chandrachur Ghosh, PhD

Department of Biosciences and Bioengineering

Indian Institute of Technology Roorkee

Roorkee 247667, Uttarakhand, India
cghosh@bt.iitr.ac.in, 94.c.gsh@gmail.com
