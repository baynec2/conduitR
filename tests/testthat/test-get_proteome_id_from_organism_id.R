test_that("this works with problematic ids", {
    problematic_id = 1218087
    out = expect_no_error(get_proteome_id_from_organism_id(problematic_id))
    })

