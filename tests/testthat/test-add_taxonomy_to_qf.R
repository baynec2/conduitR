# Build a small synthetic QFeatures with a peptides assay (Protein.Ids,
# Protein.Group on rowData) and a matching protein_groups assay, then run
# add_taxonomy_to_qf and verify peptide LCAs, protein-group labels, and
# rank assay column sums.

make_taxonomy_test_qf <- function() {
  # Peptides: rownames = stripped sequence; rowData has Protein.Ids + Protein.Group.
  #
  # PG_Ecoli  (one group, two peptides, both E. coli specific)
  #   pep_eco_a: Protein.Ids = P_Ec1;P_Ec2          -> species E. coli
  #   pep_eco_b: Protein.Ids = P_Ec1               -> species E. coli
  #
  # PG_EnteroMix (one group, two peptides, mixed taxonomic resolution)
  #   pep_mix_specific:  P_Ec1                     -> species E. coli
  #   pep_mix_family:    P_Ec1;P_Sa1;P_Kp1         -> family Enterobacteriaceae,
  #                                                   genus/species = NA
  #
  # PG_CrossDomain (one group, one peptide spanning two domains)
  #   pep_cross: P_Ec1;P_Hum1                      -> all NA
  pep_rows <- tibble::tribble(
    ~Stripped.Sequence,    ~Protein.Ids,         ~Protein.Group,
    "pep_eco_a",           "P_Ec1;P_Ec2",        "PG_Ecoli",
    "pep_eco_b",           "P_Ec1",              "PG_Ecoli",
    "pep_mix_specific",    "P_Ec1",              "PG_EnteroMix",
    "pep_mix_family",      "P_Ec1;P_Sa1;P_Kp1",  "PG_EnteroMix",
    "pep_cross",           "P_Ec1;P_Hum1",       "PG_CrossDomain"
  )
  n_samp <- 3
  samp <- paste0("S", seq_len(n_samp))
  pep_mat <- matrix(
    c(10, 20, 30,    # pep_eco_a
      40, 50, 60,    # pep_eco_b
       7,  8,  9,    # pep_mix_specific
      70, 80, 90,    # pep_mix_family
       1,  2,  3),   # pep_cross
    nrow = nrow(pep_rows), byrow = TRUE,
    dimnames = list(pep_rows$Stripped.Sequence, samp)
  )
  pep_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = pep_mat),
    rowData = S4Vectors::DataFrame(
      row.names     = pep_rows$Stripped.Sequence,
      Protein.Ids   = pep_rows$Protein.Ids,
      Protein.Group = pep_rows$Protein.Group
    )
  )

  # Protein groups: one row per unique Protein.Group, with nominal intensities.
  pg_names <- unique(pep_rows$Protein.Group)
  pg_mat <- matrix(
    seq_len(length(pg_names) * n_samp) * 100,
    nrow = length(pg_names), ncol = n_samp,
    dimnames = list(pg_names, samp)
  )
  pg_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = pg_mat),
    rowData = S4Vectors::DataFrame(row.names = pg_names)
  )

  QFeatures::QFeatures(list(peptides = pep_se, protein_groups = pg_se))
}

make_taxonomy_test_annotation <- function() {
  tibble::tribble(
    ~protein_id, ~domain,     ~kingdom,        ~phylum,         ~class,                ~order,             ~family,              ~genus,        ~species,
    "P_Ec1",     "Bacteria",  NA_character_,   "Proteobacteria","Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia coli",
    "P_Ec2",     "Bacteria",  NA_character_,   "Proteobacteria","Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia coli",
    "P_Sa1",     "Bacteria",  NA_character_,   "Proteobacteria","Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Salmonella",  "Salmonella enterica",
    "P_Kp1",     "Bacteria",  NA_character_,   "Proteobacteria","Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Klebsiella",  "Klebsiella pneumoniae",
    "P_Hum1",    "Eukaryota", "Metazoa",       "Chordata",      "Mammalia",            "Primates",         "Hominidae",          "Homo",        "Homo sapiens"
  )
}

test_that("per-peptide LCA reflects Protein.Ids taxonomy", {
  qf <- make_taxonomy_test_qf()
  ann <- make_taxonomy_test_annotation()
  qf <- add_taxonomy_to_qf(qf, ann)

  pep_rd <- SummarizedExperiment::rowData(qf[["peptides"]])

  # All-E.coli peptides: species-level LCA
  expect_equal(pep_rd["pep_eco_a", "species"], "Escherichia coli")
  expect_equal(pep_rd["pep_eco_a", "lca"],     "Escherichia coli")
  expect_equal(pep_rd["pep_eco_b", "species"], "Escherichia coli")
  expect_equal(pep_rd["pep_eco_b", "lca"],     "Escherichia coli")

  # Multi-genus same-family peptide: family-level LCA, genus/species NA
  expect_true(is.na(pep_rd["pep_mix_family", "species"]))
  expect_true(is.na(pep_rd["pep_mix_family", "genus"]))
  expect_equal(pep_rd["pep_mix_family", "family"], "Enterobacteriaceae")
  expect_equal(pep_rd["pep_mix_family", "lca"],    "Enterobacteriaceae")

  # Cross-domain peptide: all NA
  for (rk in c("domain","kingdom","phylum","class","order","family","genus","species","lca")) {
    expect_true(is.na(pep_rd["pep_cross", rk]),
                info = paste("rank", rk, "should be NA for cross-domain peptide"))
  }
})

test_that("per-protein-group LCA is derived from peptide LCAs (na.omit consensus)", {
  qf <- make_taxonomy_test_qf()
  ann <- make_taxonomy_test_annotation()
  qf <- add_taxonomy_to_qf(qf, ann)

  pg_rd <- SummarizedExperiment::rowData(qf[["protein_groups"]])

  # PG_Ecoli: both peptides agree at species level
  expect_equal(pg_rd["PG_Ecoli", "species"], "Escherichia coli")
  expect_equal(pg_rd["PG_Ecoli", "lca"],     "Escherichia coli")

  # PG_EnteroMix: pep_mix_specific is E. coli at species; pep_mix_family is NA at
  # species/genus and Enterobacteriaceae at family. With na.omit consensus,
  # the group inherits the most specific positive claim from its peptides:
  # species = E. coli, genus = Escherichia.
  expect_equal(pg_rd["PG_EnteroMix", "species"], "Escherichia coli")
  expect_equal(pg_rd["PG_EnteroMix", "genus"],   "Escherichia")
  expect_equal(pg_rd["PG_EnteroMix", "family"],  "Enterobacteriaceae")
  expect_equal(pg_rd["PG_EnteroMix", "lca"],     "Escherichia coli")

  # PG_CrossDomain: only peptide is all-NA; group label is all NA
  expect_true(is.na(pg_rd["PG_CrossDomain", "species"]))
  expect_true(is.na(pg_rd["PG_CrossDomain", "lca"]))
})

test_that("per-protein-group LCA collapses to NA when peptides positively disagree", {
  # Set up a group whose two peptides claim different species at the species
  # level — should collapse to NA at species (and genus), keep family.
  qf <- make_taxonomy_test_qf()
  pep_rd <- SummarizedExperiment::rowData(qf[["peptides"]])
  # Replace pep_eco_b's Protein.Group/Protein.Ids to put it in PG_Ecoli with
  # a Salmonella-only Protein.Ids (so the group has E. coli + Salmonella claims)
  pep_rd["pep_eco_b", "Protein.Ids"] <- "P_Sa1"
  SummarizedExperiment::rowData(qf[["peptides"]]) <- pep_rd

  ann <- make_taxonomy_test_annotation()
  qf <- add_taxonomy_to_qf(qf, ann)

  pg_rd <- SummarizedExperiment::rowData(qf[["protein_groups"]])
  # Peptides positively disagree at species and genus -> NA there;
  # they agree at family (Enterobacteriaceae) and above.
  expect_true(is.na(pg_rd["PG_Ecoli", "species"]))
  expect_true(is.na(pg_rd["PG_Ecoli", "genus"]))
  expect_equal(pg_rd["PG_Ecoli", "family"], "Enterobacteriaceae")
  expect_equal(pg_rd["PG_Ecoli", "lca"],    "Enterobacteriaceae")
})

test_that("rank assays aggregate peptide intensities by per-peptide LCA", {
  qf <- make_taxonomy_test_qf()
  ann <- make_taxonomy_test_annotation()
  qf <- add_taxonomy_to_qf(qf, ann)

  # species assay: only peptides with non-NA species contribute. That's
  # pep_eco_a (10,20,30), pep_eco_b (40,50,60), pep_mix_specific (7,8,9) -> all E. coli.
  # pep_mix_family and pep_cross have NA species and drop out.
  sp <- SummarizedExperiment::assay(qf[["species"]])
  expect_true("Escherichia coli" %in% rownames(sp))
  expect_equal(unname(sp["Escherichia coli", ]), c(10+40+7, 20+50+8, 30+60+9))

  # family assay: includes pep_mix_family on top of the species peptides
  fam <- SummarizedExperiment::assay(qf[["family"]])
  expect_true("Enterobacteriaceae" %in% rownames(fam))
  expect_equal(unname(fam["Enterobacteriaceae", ]),
               c(10+40+7+70, 20+50+8+80, 30+60+9+90))

  # domain assay: cross-domain peptide still drops (its domain is NA)
  dom <- SummarizedExperiment::assay(qf[["domain"]])
  expect_true("Bacteria" %in% rownames(dom))
  expect_equal(unname(dom["Bacteria", ]),
               c(10+40+7+70, 20+50+8+80, 30+60+9+90))
  expect_false("Eukaryota" %in% rownames(dom))
})

test_that("ranks where every peptide is NA are skipped (no assay created)", {
  # kingdom is NA for all bacterial peptides in our annotation, and the
  # Eukaryota peptide collapses to NA. So kingdom should be skipped.
  qf <- make_taxonomy_test_qf()
  ann <- make_taxonomy_test_annotation()
  qf <- add_taxonomy_to_qf(qf, ann)

  expect_false("kingdom" %in% names(qf))
  # Other ranks were created
  for (rk in c("domain","phylum","class","order","family","genus","species")) {
    expect_true(rk %in% names(qf), info = rk)
  }
})

test_that("missing required columns produce a clear error", {
  qf <- make_taxonomy_test_qf()
  ann <- make_taxonomy_test_annotation()

  # Drop Protein.Ids from peptides rowData
  pep_rd <- SummarizedExperiment::rowData(qf[["peptides"]])
  pep_rd$Protein.Ids <- NULL
  SummarizedExperiment::rowData(qf[["peptides"]]) <- pep_rd

  expect_error(add_taxonomy_to_qf(qf, ann), "Protein.Ids")
})
