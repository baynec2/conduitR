test_that("818 expected results", {

  # Create a temporary directory for test output
  tmp_dir <- withr::local_tempdir()

  # Run your function — assuming it has an argument for the output directory
  expect_no_error(
    download_fasta_from_organism_ids(818,
                                     fasta_destination_fp = paste0(tmp_dir,"test.fasta"),
                                     proteome_id_destination_fp = paste0(tmp_dir,"test.txt")
                                     )

  )

  downloads = readr::read_delim(paste0(tmp_dir,"test.txt")) |>
    dplyr::select(proteome_id,organism,organism_id,download_info)

  # Right now this is a other proteome
  print(downloads)

})


test_that("9606 expected results", {

  # Create a temporary directory for test output
  tmp_dir <- withr::local_tempdir()

  # Run your function — assuming it has an argument for the output directory
  expect_no_error(
    download_fasta_from_organism_ids(9606,
                                     fasta_destination_fp = paste0(tmp_dir,"test.fasta"),
                                     proteome_id_destination_fp = paste0(tmp_dir,"test.txt")
    )

  )

  downloads = readr::read_delim(paste0(tmp_dir,"test.txt")) |>
    dplyr::select(proteome_id,organism,organism_id,download_info)

  # This should provide a reference proteome
  print(downloads)

})



test_that("83333 expected results: reference level strain", {

  # Create a temporary directory for test output
  tmp_dir <- withr::local_tempdir()

  # Run your function — assuming it has an argument for the output directory
  expect_no_error(
    download_fasta_from_organism_ids(83333,
                                     fasta_destination_fp = paste0(tmp_dir,"test.fasta"),
                                     proteome_id_destination_fp = paste0(tmp_dir,"test.txt")
    )

  )

  downloads = readr::read_delim(paste0(tmp_dir,"test.txt")) |>
    dplyr::select(proteome_id,organism,organism_id,download_info)

  # Right now this should return the reference proteome
  print(downloads)

})

test_that("8888 expected results: not valid taxa id", {

  # Create a temporary directory for test output
  tmp_dir <- withr::local_tempdir()

  # Run your function — assuming it has an argument for the output directory
  expect_no_error(
    download_fasta_from_organism_ids(8888,
                                     fasta_destination_fp = paste0(tmp_dir,"test.fasta"),
                                     proteome_id_destination_fp = paste0(tmp_dir,"test.txt")
    )

  )

  downloads = readr::read_delim(paste0(tmp_dir,"test.txt")) |>
    dplyr::select(proteome_id,organism,organism_id,download_info)

  # Right now this should return the reference proteome
  print(downloads)

})


test_that("1682 expected results: Reference Proteome", {

  # Create a temporary directory for test output
  tmp_dir <- withr::local_tempdir()

  # Run your function — assuming it has an argument for the output directory
  expect_no_error(
    download_fasta_from_organism_ids(1682,
                                     fasta_destination_fp = paste0(tmp_dir,"test.fasta"),
                                     proteome_id_destination_fp = paste0(tmp_dir,"test.txt")
    )

  )

  downloads = readr::read_delim(paste0(tmp_dir,"test.txt")) |>
    dplyr::select(proteome_id,organism,organism_id,download_info)

  # Right now this should return the reference proteome
  print(downloads)

})


test_that("works with all", {

  all_ids = c(1682,8888,83333,9606,818)

  # Create a temporary directory for test output
  tmp_dir <- withr::local_tempdir()

  # Run your function — assuming it has an argument for the output directory
  expect_no_error(
    download_fasta_from_organism_ids(all_ids,
                                     fasta_destination_fp = paste0(tmp_dir,"test.fasta"),
                                     proteome_id_destination_fp = paste0(tmp_dir,"test.txt")
    )

  )

  downloads = readr::read_delim(paste0(tmp_dir,"test.txt")) |>
    dplyr::select(proteome_id,organism,organism_id,download_info)

  # Right now this should return the reference proteome
  print(downloads)

})


