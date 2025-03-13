# Testing with one ID
test_that("get_uniprot_data works as expected ", {
  expect_equal(get_uniprot_data("P01308")$Entry, "P01308")
})
# Testing with multiple IDs
test_that("get_uniprot_data works with multiple IDs", {
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
  annotated <- get_uniprot_data(uniprot_ids_150, batch_size = 150)
  # They are returned in a different order than they originate from the API
  expect_equal(sort(uniprot_ids_150), sort(annotated$Entry))
})

# Here is a different set of IDs where it doesn't work right
test_that("this works with problematic set of ids", {
  # Was having problems with this set:
  problematic_uniprot <- c(
    "A9KMM6", "A9KMM8", "A9KMN8", "A9KMP3", "A9KMP5", "A9KMP8", "A9KMQ8",
    "A9KMS1", "A9KMT1", "A9KMT2", "A9KMT3", "A9KMT4", "A9KMT7", "A9KMU2",
    "A9KMU8", "A9KMV0", "A9KMV1", "A9KMV9", "A9KMW6", "A9KMW9", "A9KMX1",
    "A9KMX3", "A9KMX4", "A9KMX5", "A9KMY9", "A9KN09", "A9KN26", "A9KN45",
    "A9KN46", "A9KN48", "A9KN73", "A9KN89", "A9KN91", "A9KN92", "A9KN98",
    "A9KNA6", "A9KNB4", "A9KNC3", "A9KNC4", "A9KND1", "A9KND4", "A9KND8",
    "A9KND9", "A9KNE0", "A9KNE1", "A9KNE2", "A9KNE3", "A9KNE7", "A9KNE8",
    "A9KNF7", "A9KNF8", "A9KNG0", "A9KNG2", "A9KNG4", "A9KNG5", "A9KNG9",
    "A9KNI6", "A9KNI7", "A9KNJ5", "A9KNJ9", "A9KNK0", "A9KNK7", "A9KNL0",
    "A9KNL8", "A9KNL9", "A9KNM2", "A9KNM3", "A9KNM9", "A9KNN2", "A9KNN5",
    "A9KNN7", "A9KNQ1", "A9KNR4", "A9KNR5", "A9KNS7", "A9KNU4", "A9KNV1",
    "A9KNV2", "A9KNV9", "A9KNW1", "A9KNW3", "A9KNW4", "A9KNW7", "A9KNW8",
    "A9KNX3", "A9KNX4", "A9KNX6", "A9KNX7", "A9KNX9", "A9KNZ2", "A9KNZ3",
    "A9KP10", "A9KP16", "A9KP17", "A9KP19", "A9KP20", "A9KP23", "A9KP24",
    "A9KP28", "A9KP29", "A9KP38", "A9KP44", "A9KP45", "A9KP47", "A9KP50",
    "A9KP52", "A9KP66", "A9KP79", "A9KP82", "A9KP84", "A9KP85", "A9KP87",
    "A9KP97", "A9KP98", "A9KP99", "A9KPA7", "A9KPA9", "A9KPE2", "A9KPF0",
    "A9KPF1", "A9KPF2", "A9KPF3", "A9KPF5", "A9KPG1", "A9KPG2", "A9KPG3",
    "A9KPH4", "A9KPH7", "A9KPI4", "A9KPI6", "A9KPI7", "A9KPJ8", "A9KPK8",
    "A9KPL1", "A9KPL9", "A9KPP1", "A9KPP2", "A9KPP4", "A9KPP8", "A9KPQ1",
    "A9KPQ3", "A9KPQ5", "A9KPQ6", "A9KPQ9", "A9KPT4", "A9KPT5", "A9KPU1",
    "A9KPU2", "A9KPV3", "A9KPV4"
  )

  annotated <- get_uniprot_data(problematic_uniprot,
    batch_size = 150,
    columns = "accession"
  )

  # They are returned in a different order than they originate from the API
  expect_equal(sort(problematic_uniprot), sort(annotated$Entry))
})


# Test that this works with the subset of the 5000 that was causing problems
# previously.
test_that("get_uniprot_data works with a set of 5000 ids", {
  uniprot_ids_5000 <- readRDS("tests/testthat/5000_uniprot_ids.rds")
  batch_not_working <- uniprot_ids_5000[301:450]
  results <- get_uniprot_data(batch_not_working, columns = "accession")

  expect_equal(sort(batch_not_working), sort(results$Entry))
})
