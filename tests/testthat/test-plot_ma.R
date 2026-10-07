test_that("plot_ma() returns a ggplot with default significance coloring", {
  limma_stats <- tibble::tibble(
    AveExpr   = runif(50, 20, 30),
    logFC     = rnorm(50),
    adj.P.Val = runif(50)
  )
  p <- plot_ma(limma_stats)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ma() works with color_by", {
  limma_stats <- tibble::tibble(
    AveExpr   = runif(50, 20, 30),
    logFC     = rnorm(50),
    adj.P.Val = runif(50),
    group     = rep(c("A", "B"), 25)
  )
  p <- plot_ma(limma_stats, color_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ma() works with add_loess = FALSE", {
  limma_stats <- tibble::tibble(
    AveExpr   = runif(50, 20, 30),
    logFC     = rnorm(50),
    adj.P.Val = runif(50)
  )
  p <- plot_ma(limma_stats, add_loess = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ma() works with facet_formula", {
  limma_stats <- tibble::tibble(
    AveExpr   = runif(50, 20, 30),
    logFC     = rnorm(50),
    adj.P.Val = runif(50),
    batch     = rep(c("B1", "B2"), 25)
  )
  p <- plot_ma(limma_stats, facet_formula = ~batch)
  expect_s3_class(p, "ggplot")
})
