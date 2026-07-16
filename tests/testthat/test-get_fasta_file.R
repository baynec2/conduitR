test_that("this works", {
  skip_if_offline()
  tmp_dir <- tempdir()  # Automatically deleted after session
  proteome_id = "UP000270929"
  expect_no_error(get_fasta_file(proteome_id))
  unlink(tmp_dir)
})

# UP000050895 has no sequences in UniProtKB; we fall back to UniParc (may or may not have sequences).
test_that("proteome with empty UniProtKB returns valid source (uniparc or not_downloaded)", {
  skip_if_offline()
  tmp_dir <- tempdir()
  proteome_id <- "UP000050895"
  result <- get_fasta_file(proteome_id, tmp_dir)
  expect_true(result$source %in% c("uniparc", "not_downloaded"))
  unlink(tmp_dir)
})

# UniParc /search deflines are bare (">UPI... status=active"); get_fasta_file
# stamps OS=/OX= from the proteome record so extract_fasta_info() can recover the
# taxon, without disturbing the first-token protein_id that DIA-NN keys on.
test_that("UniParc fallback deflines are stamped with OX= so taxonomy resolves", {
  skip_if_offline()
  tmp_dir <- tempdir()
  proteome_id <- "UP000050895"
  result <- get_fasta_file(proteome_id, tmp_dir)
  skip_if_not(result$source == "uniparc", "proteome did not fall back to UniParc")

  fp <- file.path(tmp_dir, paste0(proteome_id, ".fasta"))
  hdr <- grep("^>", readLines(fp), value = TRUE)
  # every UniParc defline now carries OX=, and the leading token is still the UPI id
  expect_true(all(grepl("OX=", hdr)))
  expect_true(all(grepl("^>UPI\\S* ", hdr)))

  info <- extract_fasta_info(fp)
  expect_false(any(is.na(info$organism_id)))       # taxon recovered for every protein
  expect_true(all(grepl("^UPI", info$protein_id))) # protein_id unchanged (matches DIA-NN)
  unlink(fp)
})

test_that("this works with a problematic id",{
  skip_if_offline()
  # Create a temporary directory
  tmp_dir <- tempdir()  # Automatically deleted after session
  proteome_id = "UP000041370"
  # Define temporary file path
  expect_no_error(get_fasta_file(proteome_id,file.path(tmp_dir)))
  unlink(tmp_dir)
})
