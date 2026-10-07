test_that("this works with simple test set.", {
  skip_if_offline()
  # Human and B.theta
  human_btheta = c(9606,818)
  # getting the taxonomy
  out = get_ncbi_taxonomy(human_btheta)
  # Should have Human and b.theta.
  expect_contains(out$species,c("Homo sapiens","Bacteroides thetaiotaomicron"))
})
test_that("there are no NAs",{
  skip_if_offline()
  # Human and B.theta
  human_btheta <- c(9606, 818)

  # Getting the taxonomy
  out <- get_ncbi_taxonomy(human_btheta)

  # Ensure no missing values
  expect_false(any(is.na(out$domain)))
  expect_false(any(is.na(out$kingdom)))
  expect_false(any(is.na(out$phylum)))
  expect_false(any(is.na(out$class)))
  expect_false(any(is.na(out$order)))
  expect_false(any(is.na(out$family)))
  expect_false(any(is.na(out$genus)))
  expect_false(any(is.na(out$species)))

})

test_that("problematic IDs work",{
  skip_if_offline()
  problematic_ids = c(77133,2320102)

  out = get_ncbi_taxonomy(problematic_ids)
})

# --- species/strain separation (conduitR #21) ---------------------------------
# A strain taxid must contribute a *species*-rank name to `species`, with the
# infraspecific part carried in `strain`. Otherwise two strains of one species
# are treated as two species and peptides shared between them fall to genus.

test_that("strain taxids resolve to their species-rank name", {
  skip_if_offline()
  # 349741 = Akkermansia muciniphila ATCC BAA-835 (rank: strain)
  #          parent species 239935 = Akkermansia muciniphila
  # 272559 = Bacteroides fragilis NCTC 9343 (rank: strain)
  #          parent species 817 = Bacteroides fragilis
  out <- get_ncbi_taxonomy(c(349741, 272559))

  expect_equal(
    out$species[out$organism_id == "349741"],
    "Akkermansia muciniphila"
  )
  expect_equal(
    out$strain[out$organism_id == "349741"],
    "ATCC BAA-835"
  )
  expect_equal(
    out$species[out$organism_id == "272559"],
    "Bacteroides fragilis"
  )
  expect_equal(
    out$strain[out$organism_id == "272559"],
    "NCTC 9343"
  )
})

test_that("a strain and its parent species share one species name", {
  skip_if_offline()
  # The whole point of #21: these two proteomes must agree at species rank so
  # that a peptide shared between them resolves to species, not genus.
  out <- get_ncbi_taxonomy(c(349741, 239935))

  expect_length(unique(out$species), 1L)
  expect_equal(unique(out$species), "Akkermansia muciniphila")

  # Both remain individually addressable.
  expect_setequal(out$organism_id, c("349741", "239935"))
  expect_true(is.na(out$strain[out$organism_id == "239935"]))
  expect_equal(out$strain[out$organism_id == "349741"], "ATCC BAA-835")
})

test_that("species-rank taxa are unchanged and carry no strain", {
  skip_if_offline()
  out <- get_ncbi_taxonomy(c(239935, 817))

  expect_setequal(
    out$species,
    c("Akkermansia muciniphila", "Bacteroides fragilis")
  )
  expect_true(all(is.na(out$strain)))
})

test_that("informal designations are never merged", {
  skip_if_offline()
  # 397290/1898203 are both rank `species` with no species-rank ancestor, so a
  # rank-keyed rule leaves them distinct. A name-prefix rule would wrongly fold
  # 'Lachnospiraceae bacterium A2' into 'Lachnospiraceae bacterium', promoting
  # family-level signal to a species call.
  out <- get_ncbi_taxonomy(c(397290, 1898203))

  expect_length(unique(out$species), 2L)
  expect_contains(
    out$species,
    c("Lachnospiraceae bacterium A2", "Lachnospiraceae bacterium")
  )
  expect_true(all(is.na(out$strain)))
})

test_that("strain column is always present in the returned schema", {
  skip_if_offline()
  out <- get_ncbi_taxonomy(9606)
  expect_true("strain" %in% colnames(out))
})
