# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is conduitR?

An R package for **metaproteomics** — identifying and quantifying proteins from microbial communities. It powers the [Conduit](https://github.com/baynec2/conduit) Snakemake workflow and [Conduit-GUI](https://github.com/baynec2/conduit-GUI) Shiny app, but is usable standalone.

## Common Commands

```r
# Install dependencies (from R console)
devtools::install_deps()

# Load package during development
devtools::load_all()

# Run all tests (from R console)
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-diann_to_qfeatures.R")

# Generate/update documentation
devtools::document()

# Check the full package
devtools::check()
```

Tests use `test_path("fixtures/...")` to reference fixture files under `tests/testthat/fixtures/`, which works under both `devtools::test()` and `R CMD check`. Shared path helpers are in `tests/testthat/helper-fixtures.R`. The file `inst/extdata/diann.parquet` is a user-facing example and is referenced via `system.file("extdata/diann.parquet", package = "conduitR")`.

## Architecture

### Core Data Structures

**`conduit` S4 class** (`R/conduit.R`) — the top-level container with slots:
- `@QFeatures` — a `QFeatures` object holding all assays (precursors, peptides, protein_groups, and derived assays)
- `@metrics` — named list of tibbles (DIA-NN stats, taxonomy summaries, coverage)
- `@database` — tibble of all proteins in the reference FASTA with taxonomy IDs
- `@annotations` — long-format tibble of detected proteins enriched with GO/KEGG/EggNOG/CAZy annotations
- `@taxonomy` — tibble of detected taxa with full lineage and organism type

**`QFeatures` object** — from the Bioconductor `QFeatures` package. Each assay is a `SummarizedExperiment` where rows are features (precursors/peptides/protein groups) and columns are samples. `colData` holds sample metadata, `rowData` holds feature metadata.

### Typical Data Flow

1. **Database building**: `get_fasta_file()` / `get_fasta_files()` / `download_fasta_from_proteome_ids()` → `concatenate_fasta_files()` → `extract_fasta_info()`
2. **Import**: `diann_to_qfeatures(diann_parquet_fp)` — reads a DIA-NN parquet report, filters by Q-values, and builds a `QFeatures` with three linked assays: `precursors` → `peptides` → `protein_groups` (using DIA-NN's PG.MaxLFQ values directly)
3. **Alternative import**: `prepare_qfeature(sample_annotation_fp, vector_of_matrix_fps)` — builds `QFeatures` from pre-computed TSV matrices
4. **Transformation pipeline**: `replace_zero_with_na()` → `add_log_imputed_norm_assay(qf, assay, base, impute_method, norm_method)` — adds `{assay}_log{base}`, `{assay}_log{base}_imputed`, and `{assay}_log{base}_imputed_norm` assays
5. **Annotation**: `annotate_uniprot_ids()` (batched + parallel) → `add_annotation_to_qf()`, `add_go_to_qf()`, `add_taxonomy_to_qf()`
6. **Assembly**: `build_conduit_obj()` assembles all pieces into a `conduit` S4 object
7. **Analysis**: `perform_limma_analysis()`, `perform_ora()`, `perform_gsea()`, `predict_classification()`, `predict_regression()`
8. **Visualization**: `plot_volcano()`, `plot_heatmap()`, `plot_biplot()`, `plot_taxa_tree()`, `plot_sunburst()`, `plot_kegg_pathway()`, etc.

### External API Dependencies

Many functions call external APIs that require internet access:
- **UniProt REST API** (`https://rest.uniprot.org`) — FASTA downloads, protein annotations. Falls back from UniProtKB → UniParc for proteomes with no sequences in the primary database.
- **NCBI Entrez** (`rentrez`) — taxonomy lookups
- **KEGG REST API** (`KEGGREST`) — pathway/function data
- Tests touching these APIs use `skip_if_offline()` or are in `\dontrun{}` examples.

### Key Conventions

- **Assay naming**: derived assays follow the pattern `{source_assay}_log{base}_imputed_norm`
- **Parallelism**: `annotate_uniprot_ids()` and `get_fasta_files()` use `future`/`furrr` for parallel processing; plan resets to `sequential` after
- **Test data**: all fixtures live in `tests/testthat/fixtures/`; the canonical fixture is `fixtures/conduit.rds` (a `conduit` S4 object). `inst/extdata/diann.parquet` is a user-facing example accessed via `system.file()`.
- **File paths in tests**: always use `test_path("fixtures/...")` or the helper functions from `helper-fixtures.R` (e.g., `conduit_rds()`, `taxonomy_txt()`); never hardcode paths relative to the project root
- **`add_log_imputed_norm_assays()`** (plural) is a wrapper around `add_log_imputed_norm_assay()` that iterates over all assays — used widely in tests
