# Builds a minimal but structurally complete conduit object for testing.
# Covers all column names, assay names, and slot structure needed by the test
# suite — without depending on any fixture file.
make_minimal_conduit <- function() {
  n_samples  <- 12
  n_proteins <- 20
  samp <- paste0("S", seq_len(n_samples))

  sp_names <- c(
    "Mus musculus",
    "Bacteroides thetaiotaomicron",
    "Akkermansia muciniphila",
    "Escherichia coli"
  )

  # ---- colData (shared across all assays) --------------------------------
  col_data <- S4Vectors::DataFrame(
    row.names               = samp,
    microbiome_treatment    = factor(rep(c("none", "treated"), each = n_samples / 2)),
    immunotherapy_treatment = factor(rep(c("control", "treated"), times = n_samples / 2)),
    day                     = factor(rep(c("7", "14"), times = n_samples / 2)),
    tumor_classification    = factor(rep(c("responder", "non_responder"), each = n_samples / 2))
  )

  # ---- protein_group assay -----------------------------------------------
  set.seed(42)
  pg_mat <- matrix(
    runif(n_proteins * n_samples, 1e6, 1e9),
    nrow     = n_proteins,
    ncol     = n_samples,
    dimnames = list(paste0("PG", seq_len(n_proteins)), samp)
  )
  # ~20% NAs so imputation has something to do
  pg_mat[sample(length(pg_mat), floor(length(pg_mat) * 0.2))] <- NA

  pg_se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = pg_mat),
    rowData = S4Vectors::DataFrame(
      row.names                 = rownames(pg_mat),
      Protein.Names             = paste0("Protein_", seq_len(n_proteins)),
      Genes                     = paste0("Gene", seq_len(n_proteins)),
      First.Protein.Description = paste0("Description ", seq_len(n_proteins))
    ),
    colData = col_data
  )

  # ---- species assay -----------------------------------------------------
  sp_mat <- matrix(
    runif(length(sp_names) * n_samples, 1e6, 1e9),
    nrow     = length(sp_names),
    ncol     = n_samples,
    dimnames = list(sp_names, samp)
  )
  sp_mat[sample(length(sp_mat), 2)] <- NA

  sp_se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = sp_mat),
    rowData = S4Vectors::DataFrame(
      row.names = sp_names,
      species   = sp_names
    ),
    colData = col_data
  )

  # ---- QFeatures ---------------------------------------------------------
  qf <- QFeatures::QFeatures(
    list(protein_group = pg_se, species = sp_se),
    colData = col_data
  )

  # ---- @taxonomy ---------------------------------------------------------
  taxonomy <- tibble::tibble(
    organism_id           = c(10090L, 818L, 239935L, 562L),
    domain                = c("Eukaryota", "Bacteria", "Bacteria", "Bacteria"),
    kingdom               = c("Metazoa", NA_character_, NA_character_, NA_character_),
    phylum                = c("Chordata", "Bacteroidetes", "Verrucomicrobia", "Proteobacteria"),
    class                 = c("Mammalia", "Bacteroidia", "Verrucomicrobiae", "Gammaproteobacteria"),
    order                 = c("Rodentia", "Bacteroidales", "Verrucomicrobiales", "Enterobacterales"),
    family                = c("Muridae", "Bacteroidaceae", "Akkermansiaceae", "Enterobacteriaceae"),
    genus                 = c("Mus", "Bacteroides", "Akkermansia", "Escherichia"),
    species               = sp_names,
    reference             = c(TRUE, TRUE, TRUE, TRUE),
    downloaded_by_conduit = c(FALSE, FALSE, FALSE, FALSE)
  )

  # ---- @metrics ----------------------------------------------------------
  metrics <- list(
    protein_coverage_species = tibble::tibble(
      taxon               = sp_names,
      n_proteins_detected = c(100L, 50L, 30L, 80L),
      n_proteins_db       = c(200L, 100L,  60L, 160L),
      coverage            = c(50,   50,    50,   50)
    )
  )

  # ---- @database (minimal) -----------------------------------------------
  database <- tibble::tibble(
    protein_id  = paste0("P", seq_len(n_proteins)),
    organism_id = rep(c(10090L, 818L, 239935L, 562L), times = c(5, 5, 5, 5))
  )

  # ---- @annotations (minimal) --------------------------------------------
  annotations <- tibble::tibble(
    Protein.Group = paste0("PG", seq_len(10)),
    protein_id    = paste0("P",  seq_len(10)),
    species       = rep(sp_names, times = c(3, 3, 2, 2))
  )

  new("conduit",
      QFeatures   = qf,
      metrics     = metrics,
      database    = database,
      annotations = annotations,
      taxonomy    = taxonomy,
      provenance  = NULL)
}
