

test10 = tibble::tibble(organism_type = paste0("host",1:10),
                        domain = paste0("domain",1:10),
                        kingdom = paste0("kingdom",1:10),
                        phylum = paste0("phylum",1:10),
                        class = paste0("class",1:10),
                        order = paste0("order",1:10),
                        family = paste0("family",1:10),
                        genus = paste0("genus",1:10),
                        species = paste0("species",1:10))

test_that("plot_taxa_tree() works", {
plot_taxa_tree(test10,node_color = "domain")
})

# Was having some problems for a case where there is ony one taxa

test1 = tibble::tibble(organism_type = paste0("host",1),
                        domain = paste0("domain",1),
                        kingdom = paste0("kingdom",1),
                        phylum = paste0("phylum",1),
                        class = paste0("class",1),
                        order = paste0("order",1),
                        family = paste0("family",1),
                        genus = paste0("genus",1),
                        species = paste0("species",1))

test_that("plot_taxa_tree() works", {
  plot_taxa_tree(test1, node_size = NULL)
})

test2 = tibble::tibble(organism_type = paste0("host",1:2),
                       domain = paste0("domain",1:2),
                       kingdom = paste0("kingdom",1:2),
                       phylum = paste0("phylum",1:2),
                       class = paste0("class",1:2),
                       order = paste0("order",1:2),
                       family = paste0("family",1:2),
                       genus = paste0("genus",1:2),
                       species = paste0("species",1:2))

test_that("plot_taxa_tree() works", {
  expect_no_error(plot_taxa_tree(test2, node_size = NULL))
})


test_that("plot_taxa_tree() works with info column", {
  taxonomy_data = tibble::tibble(organism_type = paste0("host",1:2),
                       domain = paste0("domain",1:2),
                       kingdom = paste0("kingdom",1:2),
                       phylum = paste0("phylum",1:2),
                       class = paste0("class",1:2),
                       order = paste0("order",1:2),
                       family = paste0("family",1:2),
                       genus = paste0("genus",1:2),
                       species = paste0("species",1:2),
                       info = c("avalible","not_avalible"))
  expect_no_error(plot_taxa_tree(taxonomy_data,node_color = "info"))
})

test_that("this works with real data",{
taxonomy_data = readr::read_delim("tests/data/01_taxonomy.txt")
expect_no_error(plot_taxa_tree(taxonomy_data,node_color = "download_info"))
})


test_that("this works with real data",{
  taxonomy_data = readr::read_delim("tests/data/01_taxonomy.txt")
  expect_no_error(plot_taxa_tree(taxonomy_data,node_size = NULL))
})
