test_that("this works with problematic ids", {
  skip_if_offline()
    problematic_id = 1218087
    out = expect_no_error(get_proteome_id_from_organism_id(problematic_id))
    })

# Regression: the UniProt `proteome_type` query field is a string enum
# ("REFERENCE"/"NON_REFERENCE"/"EXCLUDED"), not a numeric code. An earlier
# implementation queried `proteome_type:1..4`, which silently matched NOTHING
# after UniProt changed the enum -- so EVERY organism (even ones with a clear
# reference proteome) resolved to a NA proteome_id, collapsing the whole
# database down to just the directly-appended host. This must be a LIVE test:
# only a real API call detects a future enum change. A green mock would happily
# keep asserting a stale query that no longer matches anything.
test_that("an organism with a reference proteome resolves to a non-NA proteome id (proteome_type enum regression)", {
  skip_if_offline()
  out <- get_proteome_id_from_organism_id(9606)  # Homo sapiens
  expect_false(is.na(out$proteome_id),
               info = "9606 returned NA -- suspect a UniProt proteome_type enum/API change")
  expect_match(out$proteome_id, "^UP[0-9]+$")
  expect_equal(out$proteome_id, "UP000005640")  # stable human reference proteome
})
