# Methods

## Study overview

This repository analyzes the molecular composition of the mouse anterior olfactory nucleus (AON) and evaluates a predefined marker hypothesis for a dorsolateral glutamatergic population. It combines three analytically distinct resources:

1. One in-house 10x single-nucleus RNA-seq library for quality control, clustering, lineage annotation, and within-sample candidate-signature localization;
2. Allen Brain Cell Atlas (ABC Atlas) MERFISH data for spatial localization of AON glutamatergic transcriptomic clusters; and
3. ABC Atlas 10x Whole Mouse Brain data for donor-level differential expression and independent reference concordance.

The analysis does **not** observe axonal projection status. Dorsolateral location is used as a spatial proxy motivated by prior tracing literature. Results are therefore described as a **dorsolateral molecular signature**, not as a direct transcriptomic measurement of contralaterally projecting neurons.

All repository-relative paths and analysis thresholds are defined in [`config/config.yaml`](../config/config.yaml). The core workflow is orchestrated with Snakemake. Detailed outputs are written to [`results/tables/`](../results/tables/) and [`results/figures_raw/`](../results/figures_raw/); curated panels are written to [`results/figures_final/`](../results/figures_final/).

## Data sources and analytical roles

| Resource | Scope used here | Analytical role | Replication and inference boundary |
|---|---|---|---|
| In-house AON 10x snRNA-seq | One library from one mouse AON hemisphere | QC, ambient-RNA correction, doublet removal, clustering, lineage annotation, and in-house candidate-signature localization | Descriptive single-mouse analysis; no population-level or between-condition inference |
| ABC Atlas MERFISH, specimen `C57BL6J-638850` | 27,076 cells from 20 AON glutamatergic clusters across five supertypes | Spatial definition of anteroposterior and dorsolateral axes | Spatial reference; projection status is not observed |
| ABC Atlas 10x Whole Mouse Brain, `WMB-10Xv2-OLF` | 28,330 AON glutamatergic cells after scope restriction | Paired donor-level pseudobulk differential expression | Ten paired donors are the units of replication |
| ABC Atlas 10x OLF reference | 171,284 labelled cells before subclass-size filtering | Correlation-based reference label transfer | Reference concordance, not ground-truth annotation |

## Reproducibility and software environment

### Workflow organization

The core Snakemake DAG executes:

```text
preprocess -> cluster -> abca_spatial -> de -> enrich -> concordance -> crosscheck
```

The rules execute notebooks in place with `jupyter nbconvert`. The independent differential-expression cross-check is run in R with [`workflow/scripts/deseq2_crosscheck.R`](../workflow/scripts/deseq2_crosscheck.R).

Two additional stages are outside the core DAG:

- [`notebooks/07_scanvi_label_transfer.ipynb`](../notebooks/07_scanvi_label_transfer.ipynb), an optional GPU-based scANVI analysis; and
- [`notebooks/08_figure_composition.py`](../notebooks/08_figure_composition.py), which creates the curated final figure set after the required analysis outputs exist.

### Configuration

The following parameter groups are centralized in [`config/config.yaml`](../config/config.yaml):

- input and output paths;
- nucleus- and gene-level QC thresholds;
- DecontX iteration limit;
- highly variable gene count, PCA dimensions, neighbour-graph settings, and Leiden resolution sweep;
- lineage-annotation thresholds;
- spatial boundary values;
- pseudobulk cell-count and expression filters;
- differential-expression thresholds; and
- enrichment libraries.

### Environment specification

[`environment.yml`](../environment.yml) declares Python 3.12 and version-constrained dependencies including Scanpy, AnnData, scikit-learn, Snakemake, pydeseq2, DecontX, Scrublet, scvi-tools, R, and Bioconductor DESeq2. The `abc_atlas_access` dependency is pinned to a specific Git commit. Most remaining packages use minimum-version constraints rather than exact resolved versions, so the environment file is **version-constrained rather than fully locked**.

A [`Dockerfile`](../Dockerfile) builds the Conda environment in a micromamba base image. The container is intended to provide environment parity, but it is not currently built or executed by continuous integration.

### Continuous integration

The GitHub Actions workflow runs Ruff, validates YAML files, checks that core workflow files are tracked, and performs a Snakemake dry-run with placeholder inputs. It verifies repository structure and DAG resolution, not scientific reproducibility of the full data-intensive run.

## In-house preprocessing

### Input matrix

The in-house input is the Cell Ranger `filtered_feature_bc_matrix.h5` file placed at the path configured by `paths.aon_10x`. The filtered Cell Ranger matrix is used as the starting point; empty droplets are not re-called in this repository.

### Quality-control metrics and filtering

Per-nucleus QC metrics are computed before thresholds are applied. The observed distributions are plotted and compared with a robust median plus or minus three median absolute deviations as a diagnostic reference. The final thresholds remain explicit, fixed values in the configuration file.

| Parameter | Threshold | Purpose |
|---|---:|---|
| Detected genes per nucleus | at least 500 | remove low-complexity droplets and debris |
| Detected genes per nucleus | fewer than 6,000 | remove high-complexity profiles enriched for potential multiplets |
| Mitochondrial count fraction | less than 2% | remove damaged or contaminated nuclear profiles |
| Nuclei per gene | at least 3 | remove extremely sparse genes |

The stringent mitochondrial threshold reflects the nuclear assay: intact nuclei contain little cytoplasm, so elevated mitochondrial signal is treated as a damage or contamination indicator.

### Ambient-RNA correction

Ambient RNA is estimated with DecontX after a preliminary Leiden partition of the QC-filtered raw counts. Cluster labels are supplied to the model so contamination is estimated relative to broad expression populations.

The configured run uses `max_iter = 120` and reaches the iteration limit before satisfying the convergence criterion. Consequently:

- the contamination values are not interpreted as converged quantitative estimates;
- no claim is made that the final values are formal upper or lower bounds;
- raw counts are retained in `layers['counts_raw']`;
- the DecontX-adjusted matrix is stored in `.X` and `layers['counts']`; and
- downstream lineage calls and marker fractions are evaluated in a DecontX sensitivity analysis.

The post-DecontX checkpoint contains 35,028 nuclei and is loaded by the clustering notebook.

## Doublet detection, normalization, and dimensionality reduction

### Doublet scoring

Scrublet is run on the retained **raw counts**, not on the DecontX-adjusted matrix. Decontamination can shrink expression profiles toward cluster centroids and reduce the variation used by doublet simulation and scoring.

The automatic Scrublet threshold falls in the extreme score tail and does not separate a clear bimodal distribution. The analysis therefore uses the score cutoff corresponding to the assay's expected approximately 6% multiplet rate. This removes 2,076 nuclei and leaves 32,952 nuclei for clustering.

### Normalization and feature selection

Counts are normalized to the median library size and transformed with `log1p`. The 2,000 most highly variable genes are selected for dimensionality reduction.

PCA is computed with 50 components, of which the first 30 are retained for graph construction and clustering.

### Neighbour graphs and UMAP

Two graphs are constructed in the same 30-PC space:

- a 15-neighbour graph used by Leiden clustering; and
- a separate 50-neighbour graph used only for UMAP visualization, with `min_dist = 0.8`.

Separating the graphs allows a more connected, legible embedding without changing the graph partition used for cluster inference.

## Leiden resolution selection

Leiden clustering is evaluated across resolutions `0.1, 0.2, 0.3, 0.5, 0.8, 1.0`. Each solution is scored by silhouette width in PCA space using a subsample for computational tractability. Solutions are also required to maintain a minimum cluster size of 20 nuclei.

Resolution 0.1 is selected and recorded in the configuration file. It produces 22 clusters and lies on a stable silhouette plateau spanning resolutions 0.1 through 0.3. Post-clustering diagnostics test whether any retained cluster is primarily explained by doublet score or estimated contamination.

## Lineage annotation

Clusters are annotated from the **fraction of nuclei expressing** canonical marker panels rather than from mean expression alone. Fraction-expressing summaries reduce sensitivity to a small number of high-expression nuclei.

For each cluster:

1. the lineage with the largest marker-positive fraction is selected;
2. the dominant lineage must be expressed in at least 25% of nuclei;
3. a second lineage is reported as co-dominant only when it also exceeds 25% and lies within 15 percentage points of the dominant lineage; and
4. clusters that clear no lineage threshold remain unassigned.

This procedure identifies 8 excitatory clusters containing 8,764 nuclei and 13 inhibitory clusters containing 23,851 nuclei. One cluster containing 337 nuclei remains unassigned and is retained for follow-up.

The neuronal excitatory or inhibitory call is cross-checked with a mean-log-normalized-expression argmax. The two summaries agree across the called neuronal clusters. The large inhibitory fraction is supported by Gad-family expression and by reference matches to olfactory-tubercle, striatal, interneuron, and basal-forebrain classes. Its exact proportion is still interpreted in the context of the dissection boundary, which includes tissue adjacent to the AON.

## In-house candidate-signature analysis

The predefined eight-gene candidate set is:

```text
Robo2, Abi3bp, Gabrg1, Adcyap1, Chrm3, Rprm, Thrb, Cntn5
```

Expression for each candidate is standardized and averaged into a composite eight-gene score. Cluster-level summaries are used to identify the excitatory population most similar to the hypothesized dorsolateral population.

The strongest in-house score localizes to excitatory cluster 6. Because the in-house dataset contains one mouse and no retrograde projection label, this result is described only as candidate-signature localization within the sample. It is not interpreted as proof of projection identity.

## ABC Atlas MERFISH spatial analysis

The in-house snRNA-seq data contain no spatial coordinates. Spatial organization is therefore assessed with ABC Atlas MERFISH data for five AON glutamatergic supertypes, `IT AON-TT-DP Glut_1` through `IT AON-TT-DP Glut_5`.

Two axes are defined independently and reported according to their distinct interpretation.

### Anteroposterior axis

The Atlas supertypes are transcriptomically defined independently of their MERFISH coordinates. Their relationship to anteroposterior position can therefore be tested without defining groups from the same coordinate being tested.

A Kruskal-Wallis test compares cluster-level Z-position distributions. Effect size is reported with epsilon-squared, and classification separation is summarized with AUROC.

The posterior or anterior supertype boundary is the midpoint between the median Z positions of the two supertypes that straddle the boundary. With the configured threshold of 11.6:

- `Glut_3` and `Glut_5` are labelled posterior; and
- `Glut_1`, `Glut_2`, and `Glut_4` are labelled anterior.

The assignment is exported to [`results/tables/spatial_group_assignment.csv`](../results/tables/spatial_group_assignment.csv).

The anteroposterior axis is the stronger molecular axis in these data, with epsilon-squared 0.73 and AUROC 0.96.

### Dorsolateral versus ventromedial axis

The projection-motivated axis is based on a coordinate-derived dorsolateral score. For each MERFISH cell:

1. lateral distance is calculated as the absolute X displacement from the median X coordinate within the same section;
2. dorsal position is represented as negative Y;
3. the lateral and dorsal terms are standardized separately; and
4. the standardized terms are summed.

Formally,

```text
lateral_distance = |x - median_section(x)|
dorsal_term      = -y
dorsolateral_score = z(lateral_distance) + z(dorsal_term)
```

The bilateral section-specific midline prevents left- and right-hemisphere cells from cancelling one another. Clusters are ranked by median dorsolateral score, and the upper half of the 20 clusters is labelled dorsolateral. The mapping is exported to [`results/tables/spatial_cluster_dl_assignment.csv`](../results/tables/spatial_cluster_dl_assignment.csv).

The dorsolateral axis has epsilon-squared 0.32 and AUROC 0.80. It is correlated with, but not equivalent to, the anteroposterior axis. Approximately 17% of cells receive discordant group labels across the two schemes.

The coordinate split is treated as a reproducible spatial proxy. Projection identity is not observed, and cortical soma position does not necessarily determine projection target. A definitive test would require retrograde-labelled single-nucleus sequencing.

## ABC Atlas pseudobulk differential expression

### Data restriction and joins

The primary differential-expression analysis uses the ABC Atlas `WMB-10Xv2-OLF` raw-count matrix and the AON glutamatergic clusters assigned by the MERFISH spatial analysis.

Three design corrections are applied before testing:

1. **Unique join key:** expression rows and metadata are joined on `cell_label`, not the reused 16-base cell barcode.
2. **Biological replication:** raw counts are aggregated by donor and spatial group, so the donor is the inferential unit.
3. **Single chemistry:** the analysis is restricted to the 10x v2 OLF matrix to avoid a cross-chemistry batch effect.

### Pseudobulk construction

Raw counts are summed for each `(donor, group)` profile. Profiles containing fewer than 10 cells are excluded. Ten donors retain sufficient cells in both the dorsolateral and ventromedial groups and enter the paired analysis.

Genes are retained when summed counts reach at least 10 in both groups. This leaves 19,024 genes for the primary test.

### Statistical model

The primary pydeseq2 model is:

```text
~ donor + group
```

The donor term absorbs factors that are constant within a donor, including sex, while the group coefficient estimates the within-donor dorsolateral versus ventromedial contrast.

At this sample size, the parametric dispersion-trend fit does not converge and pydeseq2 uses its mean-based fallback trend. This behavior is reported explicitly. Differential genes are defined with both:

- Benjamini-Hochberg FDR < 0.05; and
- absolute log2 fold change > 1.

The primary model identifies 1,347 genes, including 550 higher in dorsolateral profiles and 797 higher in ventromedial profiles.

### Candidate-marker evaluation

Two predefined candidate collections are evaluated separately.

#### Set B: eight projection-hypothesis candidates

All eight effect estimates are positive in the dorsolateral direction. Six genes pass FDR < 0.05, and five pass both FDR < 0.05 and absolute log2 fold change > 1.

A Fisher exact test evaluates whether the thresholded candidates are over-represented among all tested genes. The resulting odds ratio is 21.9 with `p = 8.27e-5`.

The genes are reported with evidence tiers:

- `Adcyap1` has direct primary-literature support for AON glutamatergic expression;
- `Robo2` is an established axon-guidance receptor, but its canonical AON role concerns the ipsilateral lateral olfactory tract rather than a commissural pathway; and
- the remaining genes are treated as prior hypotheses rather than established projection markers.

#### Set A: broader 118-gene panel

The broader neuromodulation and interneuron panel contains 118 listed genes, of which 108 are present in the tested background. Twenty-seven pass both primary thresholds. Fisher over-representation yields an odds ratio of 4.4.

The candidate results are also compared with an established cortical projection-neuron panel. The mixed effect directions do not support a coherent canonical cortical projection-class program. The eight genes are therefore interpreted as a dorsolateral molecular signature rather than a direct projection readout.

## Differential-expression robustness analyses

The primary result is evaluated with several prespecified or transparent alternatives:

- **Unpaired design:** compared with the paired hit set for sign and set concordance.
- **Minimum-cell sweep:** the pseudobulk profile floor is varied across 10, 20, and 30 cells.
- **Anteroposterior contrast:** tests whether the candidate direction persists on the dominant molecular axis.
- **Prior hand-labelled grouping:** compares the coordinate-defined split with an earlier Target or Non-target assignment.
- **Independent R refit:** reconstructs the model with Bioconductor DESeq2 from exported pseudobulk counts and metadata.
- **Cell-level Wilcoxon test:** retained only as a transparent comparison and not treated as inferential evidence because cells from one donor are not independent replicates.

The minimum-cell sweep yields 1,347, 1,487, and 1,487 thresholded genes. The anteroposterior result is 90.3% sign-concordant genome-wide with the dorsolateral contrast and has a Jaccard overlap of 0.60 among shared thresholded genes. The prior hand-labelled grouping is 83.7% sign-concordant with a Jaccard overlap of 0.49. The independent R/Bioconductor refit overlaps the Python hit set at Jaccard 0.71.

Diagnostics include a mean-average plot, per-donor pseudobulk PCA, and a mean-variance trend plot.

## Functional characterization

Two supporting enrichment analyses are performed.

### Candidate-set over-representation

Fisher exact tests quantify over-representation of Set B and Set A among the thresholded dorsolateral differential genes. The tested gene universe is the set of genes entering the primary pydeseq2 model.

### Pathway enrichment

Enrichr over-representation analysis is run through gseapy using:

- GO Biological Process 2021;
- KEGG 2019 Mouse; and
- Reactome 2022.

Pathway analysis is applied to the anteroposterior contrast because that axis captures the dominant regional molecular structure. The anterior-enriched set yields neuroactive ligand-receptor, amine-receptor signalling, and cAMP-related terms. The posterior-enriched set has no terms at FDR < 0.05; this null result is reported without assigning an unsupported pathway identity.

Adjusted p-values are calculated within each Enrichr library, so ranks are interpreted within a library rather than across libraries.

## Reference concordance

### Correlation-based transfer

Each in-house cluster and each ABC Atlas OLF subclass is reduced to a mean expression profile over the 2,000 shared highly variable genes. ABC subclasses with fewer than 50 reference cells are excluded.

In-house clusters are matched to ABC subclasses by Spearman correlation. Rank correlation is used because the in-house and reference matrices have different transformations: median-normalized `log1p` values for the in-house object and log2 values for the reference.

The transferred neurotransmitter class agrees with the in-house call for 19 of 21 neuronal clusters, representing 88% of neuronal nuclei. The adjusted Rand index is 0.64. Most disagreement is concentrated in one large cluster rather than distributed across the atlas.

The major discordant cluster is retained as inhibitory because it expresses `Gad2`, `Foxp2`, and `Meis2` and lacks the defining marker of its highest-correlating glutamatergic reference subclass. This illustrates a limitation of nearest-reference correlation when the reference lacks a well-matched AON GABAergic class.

### Probabilistic scANVI transfer

The supplementary scANVI notebook transfers ABC OLF subclass labels at the cell level and records posterior confidence. Low-confidence cells can therefore be marked as unassignable rather than force-mapped.

This GPU step is outside the Snakemake DAG. Its cluster-level confidence summaries are used as a supplementary check on the correlation-based concordance analysis.

## Figure generation and accessibility

Notebook-level figures use a shared style module, [`notebooks/plotting_style.py`](../notebooks/plotting_style.py). The figure system applies:

- Okabe-Ito colors for low-cardinality categorical comparisons;
- a larger categorical palette for the 22 Leiden clusters, with cluster numbers placed on the embedding;
- perceptually uniform continuous color maps;
- direct labels where practical;
- visible axes, units, and threshold lines; and
- 300-dpi PNG export.

Curated multi-panel figures are generated by [`notebooks/08_figure_composition.py`](../notebooks/08_figure_composition.py). Raw diagnostics remain available separately so presentation edits do not replace the underlying analysis outputs.

## Principal output files

| Output | Description |
|---|---|
| `results/tables/retained_qc_metrics.csv` | QC metrics for retained in-house nuclei |
| `results/tables/cluster_marker_validity.csv` | cluster marker summaries and flags |
| `results/tables/spatial_group_assignment.csv` | anterior or posterior supertype assignments |
| `results/tables/spatial_cluster_dl_assignment.csv` | dorsolateral or ventromedial cluster assignments |
| `results/tables/pseudobulk_counts.csv` | donor-by-group pseudobulk count matrix |
| `results/tables/pseudobulk_coldata.csv` | donor, group, and sex metadata for pseudobulk profiles |
| `results/tables/de_dorsolateral_vs_ventromedial_deseq2.csv` | complete primary pydeseq2 results |
| `results/tables/de_significant_dorsolateral_vs_ventromedial.csv` | genes passing the primary thresholds |
| `results/tables/candidate_projection_de.csv` | eight-gene candidate result and evidence tiers |
| `results/tables/candidate_set_ora.csv` | candidate-set Fisher exact tests |
| `results/tables/deseq2_crosscheck_overlap.csv` | overlap with the independent R/Bioconductor refit |
| `results/tables/reference_concordance.csv` | correlation-based cluster-to-reference matches |
| `results/tables/scanvi_cluster_summary.csv` | cluster-level scANVI confidence summary |

## Data availability

ABC Atlas MERFISH metadata, 10x Whole Mouse Brain metadata, and `WMB-10Xv2-OLF` expression matrices are obtained from public release `20230630` with `abc_atlas_access` and placed at the paths configured in [`config/config.yaml`](../config/config.yaml).

The in-house `filtered_feature_bc_matrix.h5` is not committed. GEO or SRA deposition is pending. Large `.h5`, `.h5ad`, and `.rds` files are excluded from Git; derived tables and figures are committed where size permits.

## Limitations

1. **No direct projection measurement.** Dorsolateral location is a spatial proxy for the projection hypothesis. The definitive experiment would combine retrograde labelling with single-nucleus sequencing.
2. **Single in-house biological sample.** Clustering, annotation, and in-house candidate localization are based on one mouse hemisphere.
3. **Unequal strength of spatial axes.** Anteroposterior position explains more molecular variation than the projection-motivated dorsolateral axis.
4. **Predefined candidate panel.** Most of the eight genes are hypotheses with limited direct AON projection evidence.
5. **Reference mismatch.** The ABC Atlas does not provide a perfect matched class for every AON-adjacent inhibitory population, limiting correlation-based label transfer.
6. **DecontX non-convergence.** The configured DecontX run reaches its iteration limit. Its contamination estimates are not interpreted as converged quantitative values.
7. **Atlas DE scope.** Differential expression is limited to ten paired donors and one 10x chemistry.
8. **Environment locking.** The Conda environment uses version floors for most dependencies rather than an exact cross-platform lockfile.
9. **CI scope.** Continuous integration validates syntax, configuration, tracked workflow files, and DAG resolution but does not rerun the full analysis or build the Docker image.

## References

- Brunjes PC, Osterberg SK. Developmental markers expressed in neocortical layers are differentially exhibited in olfactory cortex. *PLOS ONE*. 2015. [doi:10.1371/journal.pone.0138541](https://doi.org/10.1371/journal.pone.0138541)
- Chen Y, Chen X, Baserdem B, et al. High-throughput sequencing of single neuron projections reveals spatial organization in the olfactory cortex. *Cell*. 2022. [doi:10.1016/j.cell.2022.09.038](https://doi.org/10.1016/j.cell.2022.09.038)
- Füllgrabe A, George N, Green M, et al. Guidelines for reporting single-cell RNA-seq experiments. *Nature Biotechnology*. 2020. [doi:10.1038/s41587-020-00744-z](https://doi.org/10.1038/s41587-020-00744-z)
- Guha T, Fertig EJ, Deshpande A. Generating colorblind-friendly scatter plots for single-cell data. *eLife*. 2022. [doi:10.7554/eLife.82128](https://doi.org/10.7554/eLife.82128)
- Illig KR, Eudy JD. Contralateral projections of the rat anterior olfactory nucleus. *Journal of Comparative Neurology*. 2009. [doi:10.1002/cne.21900](https://doi.org/10.1002/cne.21900)
- Klingler E, De la Rossa A, Fièvre S, et al. A translaminar genetic logic for the circuit identity of intracortically projecting neurons. *Current Biology*. 2019. [doi:10.1016/j.cub.2018.11.071](https://doi.org/10.1016/j.cub.2018.11.071)
- Rougier NP, Droettboom M, Bourne PE. Ten simple rules for better figures. *PLOS Computational Biology*. 2014. [doi:10.1371/journal.pcbi.1003833](https://doi.org/10.1371/journal.pcbi.1003833)
- Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021. [doi:10.1038/s41467-021-21038-1](https://doi.org/10.1038/s41467-021-21038-1)
- Wong B. Color blindness. *Nature Methods*. 2010. [doi:10.1038/nmeth1010-775a](https://doi.org/10.1038/nmeth1010-775a)
- Zhang L, et al. Behavioral role of PACAP signaling reflects its selective distribution in glutamatergic and GABAergic neuronal subpopulations. *eLife*. 2021. [doi:10.7554/eLife.61718](https://doi.org/10.7554/eLife.61718)
- Zeppilli S, Gurrola AO, Demetci P, et al. Single-cell genomics of the mouse olfactory cortex reveals contrasts with neocortex and ancestral signatures of cell type evolution. *Nature Neuroscience*. 2025. [doi:10.1038/s41593-025-01924-3](https://doi.org/10.1038/s41593-025-01924-3)
