# Document Description

This repository contains the code and processed data used for the analysis of *in vivo* CRISPR screening experiments described in the manuscript.  
All analyses were performed using the Jupyter notebooks listed below.

---

## Directory Structure

### counts/
Raw sgRNA count tables used as input for the MAGeCK analysis pipeline across all *in vivo* CRISPR screens.

### related_files/
Input files required for the data analysis described in the code notebooks.

### in_vivo_Perturb-seq/
Code and input files for Perturb-seq data analysis (related to **Fig. 4L–O** and **Fig. S3F–G**).

---

## Analysis Notebooks

### 1-MAGeCK.ipynb
Outputs of MAGeCK used for downstream data analysis.

### 2-QC.ipynb
Quality control analyses related to **Fig. 2D–G**, **Fig. S2A–D**, and **Fig. S4**.

### 3-name_conversion_to_new_mouse_gene_name_fdr.ipynb
Converts mouse gene names in the MAGeCK gene summary output to human gene names and determines hit genes based on FDR.

### 4-ROC-AUC.ipynb
Analyses related to **Fig. 2I** and **Fig. S2E**.

### 5-mageck_plot.ipynb
Generates plots related to **Fig. 3A**.

### 6-enrichment.ipynb
Enrichment analyses related to **Fig. 3F** and **Fig. 3H–J**.

### 7-compare_3_screens_complex.ipynb
Comparison of three CRISPR screens (related to **Fig. 5B**).

### 7.2-compare_3_screens-upset & PCA.ipynb
UpSet and PCA analyses comparing three screens (related to **Fig. 5C** and **Fig. S5**). 

### 8-venn-core_essential_vs_depmap_vs_invivo.ipynb
Venn diagram comparing  essential genes from DepMap, core essential genes (Hart et.al) *in vivo* screens (related to **Fig. 4F**).

### 9-compare_vglut2_vs_vgat.ipynb
Comparison between Vglut2 and Vgat screens (related to **Fig. 5D–F**).

### 10-compare_2M_vs_4M.ipynb
Comparison between 2-month and 4-month screens (related to **Fig. 5G–H**).

### 11-iNC&volcano.ipynb
Volcano plot and iNC analysis (related to **Fig. 4B–C**).

### 12-compare_invivo_vs_ineuron.ipynb
Comparison between *in vivo* and iNeuron screens (related to **Fig. 4B–C**).

### 13-4m_vs_18M&aging_hits_enrichment&compare_human_aging.ipynb
Analyses of aging-related hits and enrichment, including comparison with human aging datasets (related to **Fig. 7D–F**).

### 14-neuron-core-essential & enrichment.ipynb
Identification and enrichment analysis of neuron core essential genes (related to **Fig. 6**).

### 15-aging-sgRNA_dropout.ipynb
Analysis of sgRNA dropout in aging screens (related to **Fig. 7B**).

### 16-QC_comparing_replicates.ipynb
Quality control comparing biological replicates (related to **Fig. 2H**).

### 17-screen-pathway-label.ipynb
Pathway annotation for screen hits (related to **Fig. 3H**).

### 18-2M_4M-sgRNA_dropout.ipynb
sgRNA dropout analysis for the 2M vs 4M screens (related to **Fig. 3B**, **Fig. 3D**, and **Fig. 5B**).
