test_that("works with all organism types (reference, representative, invalid, strain)", {
  skip_if_offline()

  # Covers: reference proteome (9606), representative proteome (1682),
  # reference-level strain (83333), invalid taxa id (8888), other proteome (818)
  all_ids = c(1682, 8888, 83333, 9606, 818)

  tmp_dir <- withr::local_tempdir()

  expect_no_error(
    download_fasta_from_organism_ids(all_ids,
                                     fasta_destination_fp = paste0(tmp_dir, "test.fasta"),
                                     proteome_id_destination_fp = paste0(tmp_dir, "test.txt")
    )
  )

  downloads = readr::read_delim(paste0(tmp_dir, "test.txt"), show_col_types = FALSE) |>
    dplyr::select(proteome_id, organism, organism_id, download_info)

  print(downloads)
})
