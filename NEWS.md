# conduitR 0.1.0

First tagged release. This is the version conduit-ascent v0.1.1 and conduit-summit
build against.

## Bug fixes

* `get_ncbi_taxonomy()` no longer writes strain designations into the `species`
  column. A strain taxid (e.g. 349741, *Akkermansia muciniphila* ATCC BAA-835)
  now gets its species-rank ancestor in `species` and the infraspecific remainder
  in a new `strain` column. Before, two strains of one species had different
  `species` values, so peptides shared between them fell to genus and dropped
  out of species-level aggregation.

## Features

* The `conduit` S4 class, holding a run's `QFeatures`, database, annotations,
  taxonomy and provenance.
* `diann_to_qfeatures()` turns a DIA-NN parquet report into a `QFeatures` object
  with precursor, peptide and protein-group assays.
* `calc_taxon_fdr()`: picked target-decoy FDR at the taxon level (Savitski et
  al. 2015), used by conduit-ascent's peptidotyping.
* UniProt and NCBI access: proteome FASTA downloads, proteome metadata and NCBI
  taxonomy lookups.
* Quantification helpers: zero-to-NA, log2 imputation, species-level
  normalization and relative abundance.
* Statistics: limma differential analysis, ORA and GSEA, and LASSO, random forest
  and XGBoost models.
* Visualization: volcano plots, heatmaps, PCA, taxonomic heat trees and
  sunbursts, KEGG pathway figures, and Conduit color scales and themes.
