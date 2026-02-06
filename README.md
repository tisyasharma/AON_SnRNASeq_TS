# AON Single-Nucleus RNA-Seq Analysis

Single-nucleus RNA sequencing analysis of the mouse Anterior Olfactory Nucleus (AON), integrating 10x Genomics snRNA-seq data with the Allen Brain Cell Atlas to identify molecular markers for contralaterally-projecting neurons.

## Results

- **10,966 cells** across **28 clusters** (6 excitatory, 22 inhibitory)
- Identified spatial segregation of AON glutamatergic neurons into dorsolateral (target) and ventromedial (non-target) populations
- Top differentially expressed gene: **Abi3bp** (log2 fold-change +5.02)

## Setup

```bash
pip install -r requirements.txt
```

## Data

Place the 10x filtered feature-barcode matrix (`filtered_feature_bc_matrix.h5`) in `data/aon_10x/`.
Allen Brain Atlas files belong in `data/allen_brain_atlas/`.
Large binary files (`.h5`, `.h5ad`) are excluded here. 

## Analysis Workflow

1. **AON_snRNAseq_TS.ipynb** - Quality control, normalization (Pearson residuals), PCA,
   Leiden clustering, marker gene discovery, excitatory/inhibitory classification

2. **ABCA_Analysis.ipynb** - Spatial visualization of AON neurons using MERFISH coordinates,
   definition of Target (dorsolateral) vs Non-target (ventromedial) populations

3. **10x_ClusterAnalysis.ipynb** - Differential expression analysis between Target and
   Non-target clusters using Allen Brain Atlas 10x Whole Mouse Brain dataset

## Key Findings

| Cluster Type | Count | Markers |
|--------------|-------|---------|
| Excitatory | 6 | Slc17a7 (vGlut1) |
| Inhibitory | 22 | Gad1 |

Top candidate markers for dorsolateral AON neurons (contralaterally-projecting):
- Robo2, Abi3bp, Gabrg1, Adcyap1, Chrm3, Rprm, Thrb, Cntn5
