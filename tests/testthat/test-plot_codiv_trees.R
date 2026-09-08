make_trees <- function() {
  set.seed(7)
  h <- ape::rtree(6); s <- ape::rtree(40)
  hs <- data.frame(Host = rep(h$tip.label, length.out = 40),
                   Symbiont = s$tip.label, stringsAsFactors = FALSE)
  list(h = h, s = s, hs = hs)
}

test_that("plot_codiv_trees returns a ggplot for each color mode", {
  skip_if_not_installed("ggplot2")
  d <- make_trees()
  for (cb in c("none", "host")) {
    p <- plot_codiv_trees(d$h, d$s, d$hs, color_by = cb)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
})

test_that("plot_codiv_trees supports significance coloring with results", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("adephylo")
  d <- make_trees()
  res <- codiv(d$h, d$s, d$hs, min_hosts = 4, min_symbiont_tips = 7,
               span_fraction = 0.6, permutations = 99, methods = "hommola",
               cores = 1, seed = 1, verbose = FALSE)
  p <- plot_codiv_trees(d$h, d$s, d$hs, codiv_results = res,
                        color_by = "significance", significance_threshold = 0.5)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("significance coloring requires codiv_results", {
  skip_if_not_installed("ggplot2")
  d <- make_trees()
  expect_error(
    plot_codiv_trees(d$h, d$s, d$hs, color_by = "significance"),
    "requires"
  )
})

test_that("both cladogram and branch-length layouts build", {
  skip_if_not_installed("ggplot2")
  d <- make_trees()
  p1 <- plot_codiv_trees(d$h, d$s, d$hs, use_branch_lengths = FALSE)
  p2 <- plot_codiv_trees(d$h, d$s, d$hs, use_branch_lengths = TRUE)
  expect_silent(ggplot2::ggplot_build(p1))
  expect_silent(ggplot2::ggplot_build(p2))
})
