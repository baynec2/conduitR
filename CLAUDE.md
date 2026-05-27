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

### Aggregation & Quantification Workflow

The core analysis pattern for downstream statistics and visualization:

1. **Aggregate** features to a taxonomic or functional level:
   ```r
   qf <- aggregate_assay_by_annotation(qf, i = "peptides", fcol = "genus", include_na = "drop")
   ```
   - `i`: source assay — `"peptides"` for taxonomic aggregation, `"protein_groups"` for functional
   - `fcol`: any column registered as an aggregation target (see below)
   - `include_na`: `"drop"` excludes unannotated features, `"group"` collects them as `"Unassigned"`
   - Creates a new assay named after `fcol` (e.g., `"genus"`)

2. **Check available aggregation targets** with `aggregation_targets(qf)` — returns a named list with `kind` (`"taxonomic"` or `"functional"`) and `from` (source assay). Taxonomic targets: `domain`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`. Functional targets include: `go`, `eggnog`, `eggnog_code`, `kegg_orthology`, `kegg_map_pathway`, `cazy_class`, `cazy_family`, plus `uniprot_`-prefixed variants.

3. **Relative abundance**: `add_relative_abundance_assay(qf, assay_name = "genus")` — creates `{assay}_rel_abundance` (values sum to 100 per sample).

4. **Log-transform for statistics**: `add_log_imputed_norm_assay(qf, assay = "genus", base = 2, impute_method = "MinDet", norm_method = "none")` — creates `{assay}_log2_MinDet_none`. Zeros are replaced with NA, then log-transformed, then imputed.

5. **Tidy format**: `tidy_conduit(qf, assay_name = "genus")` — returns a long-format tibble with columns: `file` (sample ID), `rowid` (feature ID), `value`, plus all `colData` and `rowData` columns.

6. **Differential analysis**: `perform_limma_analysis(qf, assay_name, formula, contrast)` — runs limma on the specified assay. Uses QFeatures-level `colData`, so custom columns added after construction (e.g., factor conversions) are visible to the model formula.

7. **Subsetting by taxonomy for organism-specific analysis**: filter the `protein_groups` assay by `rowData` columns (e.g., `phylum == "Chordata"` for host-only) before aggregating to functional terms:
   ```r
   se <- qf[["protein_groups"]]
   qf[["protein_groups"]] <- se[which(rowData(se)$phylum == "Chordata"), ]
   qf <- aggregate_assay_by_annotation(qf, i = "protein_groups", fcol = "go", include_na = "drop")
   ```

### Annotations Slot

The `@annotations` slot is a long-format tibble with columns: `Protein.Group`, `annotation_type`, `term`, `description`. Key annotation types:
- `uniprot_go` — GO terms from UniProt (has descriptions)
- `go` — GO terms from eggNOG (descriptions are NA; use `uniprot_go` to look up descriptions for the same GO IDs)
- `eggnog` — specific COG/NOG IDs (e.g., `COG1734`)
- `eggnog_code` — single-letter COG functional categories (e.g., `J`, `K`, `M`)
- `kegg_orthology`, `kegg_map_pathway` — KEGG annotations
- `cazy_class`, `cazy_family` — CAZyme annotations
- All of the above also exist with `uniprot_` prefix (from UniProt source vs eggNOG source)

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
