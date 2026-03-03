# Simple Test
test_that("Simple Test works", {
  skip_if_offline()
  organism_ids = c("9606", "818")
  expect_no_error(get_proteome_ids_from_organism_ids(organism_ids))
})

# These IDs were briefly problematic
prob =c(10090,818,927665,1339352,1401073,1339280,457389,1287488,1339350,328812,
  357276,1574265,246787,46506,28116,821,1679444,47678,1780376,371601,162156,
  820,823,1739298,338188,674529,1796613,239935,1965650,871325,77133,2024223,
  2320102,2486468,28111,2486471,349741,449673,2030927,2489214,39491,2292240,
  702446,471870,484018,626522,537013,679189,469590,563193,657314,470145,
  226186,679935,908612,997884,679937,997887,997873,657309,665953,457387,709991,
  997881,469585,997877,585544,556259,742817,868129,888743,762982,817,667015,
  188937,295405,999418,742726,999419,742727,1262753,880074,547042,1235785,
  702447,397290,1262986,1235786,397291,397288)

test_that("it works with problematic ids", {
  skip_if_offline()
  expect_no_error(get_proteome_ids_from_organism_ids(as.character(prob)))
})

# These IDs were also problematic at some point
test_that("it works with new problematic ids",{
  skip_if_offline()
    problematic_ids = c(77133,2320102)
    out = get_proteome_ids_from_organism_ids(problematic_ids)
})

# Was having a problem where some ids were silently dropped- this happens if the
# proteome for that taxa is not in the uniprot database.
test_that("you get the same number out as you put in",{
  skip_if_offline()
  # Representative mix: organisms with reference proteomes, redundant proteomes,
  # and some that may lack a proteome entirely, covering the no-silent-drop bug.
  ids <- c(
    9606,   # Homo sapiens (reference proteome)
    818,    # Bacteroides thetaiotaomicron (gut microbe, reference proteome)
    562,    # Escherichia coli K-12 (reference proteome)
    2173,   # Methanobrevibacter smithii (gut archaeon)
    1747,   # Clostridioides difficile
    79604,  # Bifidobacterium adolescentis
    2317,   # Methanosphaera stadtmanae
    457415, # Faecalibacterium prausnitzii M21/2
    1262868,# Ruminococcaceae bacterium
    1860157 # Candidatus Methanomassiliicoccales archaeon
  )
  out = get_proteome_ids_from_organism_ids(ids)
  ids_not_in_output = setdiff(as.character(ids), as.character(out$organism_id))
  expect_equal(length(ids_not_in_output), 0)
})

# Test for duplicate ids
test_that("it works with duplicate ids",{
  skip_if_offline()
  ids = c(77133,77133)
  expect_no_error(get_proteome_ids_from_organism_ids(ids))
})


test_that("it produces NAs for organism IDs with no ID", {
  skip_if_offline()

  t = readr::read_delim(organism_ids_txt())

  organism_ids <- t$organism_id


  organism_ids <- organism_ids[!is.na(organism_ids)]

  pi = get_proteome_ids_from_organism_ids(t$organism_id)

  expect_true(any(is.na(pi$proteome_id)))
})
