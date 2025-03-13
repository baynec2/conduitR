test_that("annotate_uniprot_ids works", {
  expect_equal(annotate_uniprot_ids("P01308")$Entry, "P01308")
})

test_that("We get same ids out as we put in", {
  # 150 ids
  uniprot_ids_150 <- c(
    "A9KHN2", "A9KI94", "A9KLG0", "A9KM56", "A9KM95", "A9KMA2", "A9KMQ4",
    "A9KQ65", "A9KSH6", "A9KSR8", "A9KSR9", "A9KSS9", "A9KSW2", "A9KT32",
    "A9KTE6", "A9KH99", "A9KHD0", "A9KHD3", "A9KHD9", "A9KHE4", "A9KHI3",
    "A9KHL4", "A9KHL6", "A9KHN1", "A9KHN3", "A9KHN4", "A9KHU4", "A9KHU5",
    "A9KHW3", "A9KHX2", "A9KHX3", "A9KHY0", "A9KHY1", "A9KHY5", "A9KHY6",
    "A9KHY8", "A9KI24", "A9KI25", "A9KI30", "A9KI50", "A9KI51", "A9KI75",
    "A9KI90", "A9KI93", "A9KIA0", "A9KIA6", "A9KIB0", "A9KIG5", "A9KIN4",
    "A9KIN7", "A9KIP2", "A9KIR8", "A9KIT6", "A9KIW5", "A9KJ18", "A9KJ19",
    "A9KJ24", "A9KJ36", "A9KJ78", "A9KJG1", "A9KJG6", "A9KJH2", "A9KJH3",
    "A9KJK5", "A9KJL3", "A9KJL4", "A9KJR9", "A9KJY0", "A9KJY4", "A9KK06",
    "A9KK15", "A9KK22", "A9KK74", "A9KK75", "A9KK83", "A9KK92", "A9KK94",
    "A9KK96", "A9KK97", "A9KKE7", "A9KKE9", "A9KKG9", "A9KKH7", "A9KKR4",
    "A9KKR5", "A9KKR6", "A9KKS3", "A9KKT0", "A9KKT9", "A9KKU4", "A9KKU6",
    "A9KKW6", "A9KKY1", "A9KKY6", "A9KKY7", "A9KL08", "A9KL30", "A9KL35",
    "A9KL72", "A9KL77", "A9KL78", "A9KL97", "A9KL99", "A9KLD6", "A9KLD8",
    "A9KLG1", "A9KLJ9", "A9KLK3", "A9KLK5", "A9KLK6", "A9KLL7", "A9KLM6",
    "A9KLQ9", "A9KLR0", "A9KLR1", "A9KLR4", "A9KLR6", "A9KLR7", "A9KLS2",
    "A9KLS3", "A9KLT8", "A9KLU0", "A9KLU2", "A9KLU4", "A9KLV9", "A9KLX9",
    "A9KLZ5", "A9KM11", "A9KM25", "A9KM34", "A9KM91", "A9KMA0", "A9KMA1",
    "A9KMB6", "A9KMB8", "A9KMD8", "A9KMD9", "A9KME7", "A9KMF0", "A9KMF5",
    "A9KMI1", "A9KMN9", "A9KMP0", "A9KMQ0", "A9KMQ6", "A9KMS2", "A9KMS8",
    "A9KMU1", "A9KMU9", "A9KMV3"
  )

  annotated <- annotate_uniprot_ids(uniprot_ids_150,
    batch_size = 150,
    columns = "accession"
  )

  expect_equal(sort(annotated$accession), sort(uniprot_ids_150))
})


test_that("Does this work with a large, realistic set?", {
  # 5000 ids
  uniprot_ids_5000 <- readRDS("tests/testthat/5000_uniprot_ids.rds")
  annotated <- annotate_uniprot_ids(uniprot_ids_5000,
    batch_size = 100,
    columns = "accession"
  )

  expect_equal(sort(annotated$accession), sort(uniprot_ids_5000))
})


test_that("Does this work with a large, realistic set without parallel?", {
  # 5000 ids
  uniprot_ids_5000 <- readRDS("tests/testthat/5000_uniprot_ids.rds")
  annotated <- annotate_uniprot_ids(uniprot_ids_5000,
    batch_size = 100,
    columns = "accession",
    parallel = FALSE
  )

  expect_equal(sort(annotated$accession), sort(uniprot_ids_5000))
})



test_that("parallel is faster than non parallel", {
  # 5000 ids
  uniprot_ids_5000 <- readRDS("tests/testthat/5000_uniprot_ids.rds")

  parallel <- bench::bench_time(annotate_uniprot_ids(uniprot_ids_5000,
    batch_size = 100,
    columns = "accession",
    parallel = TRUE
  ))

  single <- bench::bench_time(annotate_uniprot_ids(uniprot_ids_5000,
    batch_size = 100,
    columns = "accession",
    parallel = FALSE
  ))

  # Parallel should be faster than non-parallel
  expect_gt(single[[2]], parallel[[2]])
})


test_that({

  annotate_uniprot_ids("Q8A1G1",columns = "xref_pfam,xref_cazy,xref_eggnog")

})
