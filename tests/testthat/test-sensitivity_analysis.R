test_that("sensitivity_analysis returns the expected structure", {
  skip_if_not_installed("adephylo")
  set.seed(42)
  h <- ape::rtree(8); s <- ape::rtree(120)
  hs <- data.frame(Host = rep(h$tip.label, length.out = 120),
                   Symbiont = s$tip.label, stringsAsFactors = FALSE)

  sens <- sensitivity_analysis(
    h, s, hs,
    span_fraction_range = c(0.3, 0.5),
    min_symbiont_tips_range = c(7, 10),
    permutations = 99, min_hosts = 4, methods = c("hommola", "paco"),
    cores = 1, seed = 1, verbose = FALSE
  )

  expect_named(sens, c("sensitivity_data", "results_by_params", "heatmap_data"))
  # one row per parameter combination (2 x 2)
  expect_equal(nrow(sens$sensitivity_data), 4)
  expect_length(sens$results_by_params, 4)
  # heatmap: combinations x methods
  expect_equal(nrow(sens$heatmap_data), 4 * 2)
  expect_setequal(sens$heatmap_data$method, c("hommola", "paco"))
  # each stored result is a codiv object
  expect_true(all(vapply(sens$results_by_params,
                         function(r) inherits(r, "codiv"), logical(1))))
  # per-method significance columns present
  expect_true(all(c("hommola_n_sig", "paco_n_sig") %in%
                    colnames(sens$sensitivity_data)))
})

test_that("sensitivity_analysis validates its ranges", {
  skip_if_not_installed("adephylo")
  h <- ape::rtree(6); s <- ape::rtree(30)
  hs <- data.frame(Host = rep(h$tip.label, 5), Symbiont = s$tip.label,
                   stringsAsFactors = FALSE)

  expect_error(
    sensitivity_analysis(h, s, hs, span_fraction_range = c(0, 1.5),
                         verbose = FALSE),
    "span_fraction_range"
  )
  expect_error(
    sensitivity_analysis(h, s, hs, min_symbiont_tips_range = c(7, 7.5),
                         verbose = FALSE),
    "min_symbiont_tips_range"
  )
})
