test_that("create_provenance() returns the correct structure", {
  prov <- create_provenance(
    workflow_version  = "1.2.3",
    uniprotkb_release = "2024_05"
  )

  expect_type(prov, "list")
  expect_named(prov, c("workflow_version", "generated_date", "uniprotkb_release", "config"))
  expect_equal(prov$workflow_version, "1.2.3")
  expect_equal(prov$uniprotkb_release, "2024_05")
  expect_s3_class(prov$generated_date, "Date")
  expect_null(prov$config)
})

test_that("create_provenance() accepts a config tibble", {
  cfg <- tibble::tibble(parameter = c("fdr", "min_peptides"), value = c("0.01", "2"))
  prov <- create_provenance(
    workflow_version  = "2.0.0",
    uniprotkb_release = "2025_01",
    config = cfg
  )

  expect_equal(nrow(prov$config), 2)
  expect_equal(names(prov$config), c("parameter", "value"))
})

test_that("create_provenance() uses the supplied generated_date", {
  d <- as.Date("2025-06-15")
  prov <- create_provenance(
    workflow_version  = "1.0.0",
    uniprotkb_release = "2025_01",
    generated_date = d
  )
  expect_equal(prov$generated_date, d)
})

test_that("create_provenance() errors on bad workflow_version", {
  expect_error(
    create_provenance(workflow_version = 123, uniprotkb_release = "2024_05"),
    "is.character"
  )
  expect_error(
    create_provenance(workflow_version = c("1.0", "2.0"), uniprotkb_release = "2024_05"),
    "length"
  )
})

test_that("create_provenance() errors on bad uniprotkb_release", {
  expect_error(
    create_provenance(workflow_version = "1.0.0", uniprotkb_release = 2024),
    "is.character"
  )
})

test_that("create_provenance() errors when config lacks required columns", {
  bad_cfg <- tibble::tibble(key = "fdr", val = "0.01")
  expect_error(
    create_provenance(
      workflow_version  = "1.0.0",
      uniprotkb_release = "2024_05",
      config = bad_cfg
    )
  )
})

# Helper: build a fresh minimal conduit object for show() tests.
# The conduit.rds fixture predates the metrics slot and cannot be used with show().
make_test_conduit <- function(provenance = NULL) {
  mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list("r1", "s1"))
  se  <- SummarizedExperiment::SummarizedExperiment(assays = list(intensity = mat))
  qf  <- QFeatures::QFeatures(
    list(test = se),
    colData = S4Vectors::DataFrame(row.names = "s1")
  )
  new("conduit",
      QFeatures   = qf,
      metrics     = list(),
      database    = tibble::tibble(),
      annotations = tibble::tibble(),
      taxonomy    = tibble::tibble(),
      provenance  = provenance)
}

test_that("show() works on a conduit object with provenance = NULL", {
  obj <- make_test_conduit(provenance = NULL)
  out <- capture.output(show(obj))
  expect_true(any(grepl("Provenance", out)))
  expect_true(any(grepl("Not Available", out)))
})

test_that("show() displays provenance details when populated", {
  prov <- create_provenance(
    workflow_version  = "1.2.3",
    uniprotkb_release = "2024_05",
    config = tibble::tibble(parameter = "fdr", value = "0.01")
  )
  obj <- make_test_conduit(provenance = prov)
  out <- capture.output(show(obj))
  expect_true(any(grepl("1\\.2\\.3", out)))
  expect_true(any(grepl("2024_05", out)))
  expect_true(any(grepl("Config entries.*1", out)))
})

test_that("get_uniprotkb_release() returns a YYYY_MM string", {
  skip_if_offline()
  rel <- get_uniprotkb_release()
  expect_type(rel, "character")
  expect_match(rel, "^\\d{4}_\\d{2}$")
})
