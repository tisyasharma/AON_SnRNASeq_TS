# AON snRNA-seq pipeline driver.
# Each stage executes a notebook (or the R cross-check) in order. Parameters live in
# config/config.yaml. Run:  snakemake --cores 4 all   (or a single target, e.g.  snakemake de)
#
# Notebooks are executed in place with nbconvert so the committed .ipynb carries the run
# outputs. Heavy stages (preprocess runs ~3 h of DecontX) dominate wall-clock.

configfile: "config/config.yaml"

RESULTS      = config["paths"]["results_dir"]
POST_DECONTX = config["paths"]["aon_post_decontx"]
PREPROCESSED = config["paths"]["aon_preprocessed"]

rule all:
    input:
        POST_DECONTX,
        PREPROCESSED,
        f"{RESULTS}/figures_raw/umap_cell_types.png",
        f"{RESULTS}/tables/reference_concordance.csv",
        f"{RESULTS}/tables/spatial_cluster_dl_assignment.csv",
        f"{RESULTS}/tables/spatial_group_assignment.csv",
        f"{RESULTS}/tables/de_dorsolateral_vs_ventromedial_deseq2.csv",
        f"{RESULTS}/tables/candidate_projection_de.csv",
        f"{RESULTS}/tables/pseudobulk_counts.csv",
        f"{RESULTS}/tables/candidate_set_ora.csv",
        f"{RESULTS}/tables/enrichment_posterior_up.csv",
        f"{RESULTS}/tables/deseq2_crosscheck_overlap.csv",


def _run_nb(nb):
    # KMP_DUPLICATE_LIB_OK avoids a macOS duplicate-OpenMP-runtime abort that otherwise kills the
    # kernel at Scrublet (multiple libomp copies get loaded transitively through numba/sklearn).
    return (
        "KMP_DUPLICATE_LIB_OK=TRUE "
        "jupyter nbconvert --to notebook --execute --inplace "
        f"--ExecutePreprocessor.timeout=-1 notebooks/{nb}.ipynb"
    )


rule preprocess:
    input:
        h5=config["paths"]["aon_10x"],
    output:
        POST_DECONTX,
    shell:
        _run_nb("01_preprocessing")


rule cluster:
    input:
        POST_DECONTX,
    output:
        h5ad=PREPROCESSED,
        umap=f"{RESULTS}/figures_raw/umap_cell_types.png",
        leiden_umap=f"{RESULTS}/figures_raw/umap_leiden_clusters.png",
        proportions=f"{RESULTS}/figures_raw/celltype_proportions.png",
        markers=f"{RESULTS}/figures_raw/marker_dotplot.png",
        signature=f"{RESULTS}/figures_raw/inhouse_candidate_signature.png",
        marker_validity=f"{RESULTS}/tables/cluster_marker_validity.csv",
    shell:
        _run_nb("02_clustering")


rule concordance:
    input:
        h5ad=PREPROCESSED,
        meta=config["paths"]["allen_10x_meta"],
        log2=config["paths"]["allen_log2_matrix"],
    output:
        table=f"{RESULTS}/tables/reference_concordance.csv",
        heatmap=f"{RESULTS}/figures_raw/reference_concordance_heatmap.png",
    shell:
        _run_nb("06_reference_concordance")


rule abca_spatial:
    input:
        meta=config["paths"]["allen_merfish_meta"],
    output:
        dl=f"{RESULTS}/tables/spatial_cluster_dl_assignment.csv",
        ap=f"{RESULTS}/tables/spatial_group_assignment.csv",
        clusters=f"{RESULTS}/figures_raw/spatial_clusters.png",
        dl_fig=f"{RESULTS}/figures_raw/spatial_dorsolateral_assignment.png",
        ap_fig=f"{RESULTS}/figures_raw/spatial_posterior_anterior.png",
    shell:
        _run_nb("03_abca_spatial")


rule de:
    input:
        matrix=config["paths"]["allen_raw_matrix"],
        meta=config["paths"]["allen_10x_meta"],
        dl=f"{RESULTS}/tables/spatial_cluster_dl_assignment.csv",
        ap=f"{RESULTS}/tables/spatial_group_assignment.csv",
    output:
        de=f"{RESULTS}/tables/de_dorsolateral_vs_ventromedial_deseq2.csv",
        ap_de=f"{RESULTS}/tables/de_posterior_vs_anterior_deseq2.csv",
        cand=f"{RESULTS}/tables/candidate_projection_de.csv",
        cand_panel=f"{RESULTS}/tables/candidate_gene_de.csv",
        sig_dl=f"{RESULTS}/tables/de_significant_dorsolateral_vs_ventromedial.csv",
        sig_ap=f"{RESULTS}/tables/de_significant_posterior_vs_anterior.csv",
        pb=f"{RESULTS}/tables/pseudobulk_counts.csv",
        cd=f"{RESULTS}/tables/pseudobulk_coldata.csv",
        effectsize=f"{RESULTS}/figures_raw/candidate_gene_effectsize.png",
        volcano_dl=f"{RESULTS}/figures_raw/volcano_dorsolateral_vs_ventromedial.png",
        volcano_ap=f"{RESULTS}/figures_raw/volcano_posterior_vs_anterior.png",
        projection=f"{RESULTS}/figures_raw/projection_panel.png",
        ma=f"{RESULTS}/figures_raw/de_diagnostics_ma.png",
        pca=f"{RESULTS}/figures_raw/de_diagnostics_pca.png",
        meanvar=f"{RESULTS}/figures_raw/de_diagnostics_meanvar.png",
    shell:
        _run_nb("04_differential_expression")


rule enrich:
    input:
        dl=f"{RESULTS}/tables/de_dorsolateral_vs_ventromedial_deseq2.csv",
        ap=f"{RESULTS}/tables/de_posterior_vs_anterior_deseq2.csv",
    output:
        ora=f"{RESULTS}/tables/candidate_set_ora.csv",
        up=f"{RESULTS}/tables/enrichment_posterior_up.csv",
        down=f"{RESULTS}/tables/enrichment_anterior_up.csv",
        dotplot=f"{RESULTS}/figures_raw/enrichment_dotplot.png",
    shell:
        _run_nb("05_enrichment")


rule crosscheck:
    input:
        pb=f"{RESULTS}/tables/pseudobulk_counts.csv",
        cd=f"{RESULTS}/tables/pseudobulk_coldata.csv",
    output:
        f"{RESULTS}/tables/deseq2_crosscheck_overlap.csv",
    shell:
        "Rscript workflow/scripts/deseq2_crosscheck.R "
        "{input.pb} {input.cd} {output}"
