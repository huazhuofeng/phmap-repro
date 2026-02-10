# Seurat v5 Robust Marker Pipeline (Snakemake)

This is a reusable end‑to‑end pipeline for robust cross‑sample marker identification (pseudobulk + edgeR + cross‑sample consistency) on **any Seurat v5 object**.

## 1) Scope and Input Requirements

The pipeline can be applied to any Seurat v5 object, provided that its `meta.data` contains the following stratification variables:
- `sample`: biological replicate unit (**required**)
- `cluster_col`: cell group label (default `final_celltype`, configurable)
- `study`: batch/study covariate (default `study`, configurable)

Supported input file formats:
- `.qs` (preferred)
- `.rds`

Configure in `config/config.yaml`:
- `input`
- `outdir`
- `sample_col` / `cluster_col` / `study_col`
- Statistical threshold parameters (`min_cells`, `fdr_max`, etc.)

## 2) Output Naming and Organization (Standardized)

All results are written to `outdir` (e.g., `results/run_001/`):

- `phase4_summary.csv`: summary status for all clusters  
  - `status`: `ok / skipped / failed`
  - `n_samples_usable`, `n_markers`, etc.
- `markers__<cluster>.csv`: final marker table per cluster
- `consistency__<cluster>.csv`: cross‑sample consistency metrics per cluster
- `full_out.rds` (and optionally `full_out.qs`): full output object archive
- `run_params.rds` + `run_params.json`: pipeline parameter archive

Logs:
- Snakemake logs: `logs/pipeline.log`

## 3) How to Run

### A. Snakemake (Recommended)

```bash
cd seuratv5_robust_markers_pipeline
snakemake --use-conda -j 4
```

### B. Run R Script Directly

```bash
Rscript scripts/run_pipeline.R \
  --input /path/to/object.qs \
  --outdir results/run_001 \
  --sample_col sample \
  --cluster_col final_celltype \
  --study_col study
```

## 4) Should a New Standard Conda Environment Be Created?

**Recommended: Yes.**

Reasons:
1. Seurat/edgeR and their dependencies are version‑sensitive.
2. Decouples from the system R environment, avoiding conflicts between projects.
3. Facilitates reproducibility (the same `envs/r-marker.yml` can be reused across machines).

This project provides:
- `envs/r-marker.yml`

If you already have a stable environment, you may skip this step; however, for reproducibility and collaboration, we recommend using the provided standard environment.

## 5) Pipeline Logic Overview

1. **Phase0**: metadata validation (missing values, consistency, sample‑to‑study mapping)
2. **Phase1**: build `(sample, cluster)` pseudobulk matrix (genes × pseudo‑samples)
3. **Phase2**: for each cluster, perform `A vs rest (within‑sample)` differential analysis using edgeR
4. **Phase3**: compute cross‑sample consistency metrics (`median_logFC_by_sample`, `prop_same_sign`, `expr_prop_in_A`) and robust ranking
5. **Phase4**: batch processing, output summary for all clusters
