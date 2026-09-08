test_that("codiv_null_scans returns a calibrated null for codiversifying data", {
  skip_if_not_installed("adephylo")
  set.seed(7)
  host <- ape::rcoal(14); host$tip.label <- paste0("H", 1:14)
  sym <- host; sym$tip.label <- paste0("S", 1:14)
  sym$edge.length <- sym$edge.length * runif(length(sym$edge.length), 0.7, 1.3)
  df <- data.frame(Host = paste0("H", 1:14), Symbiont = paste0("S", 1:14),
                   stringsAsFactors = FALSE)

  obs <- codiv(host, sym, df, min_hosts = 4, min_symbiont_tips = 7,
               span_fraction = 0.95, permutations = 99, methods = "hommola",
               cores = 1, seed = 1, subtree_features = FALSE, verbose = FALSE)

  null <- codiv_null_scans(host, sym, df, n_permutations = 15, observed = obs,
                           stat_threshold = 0.75, p_threshold = 0.05,
                           min_hosts = 4, min_symbiont_tips = 7,
                           span_fraction = 0.95, permutations = 99,
                           methods = "hommola", cores = 1, seed = 1,
                           verbose = FALSE)

  expect_named(null, c("observed", "null_counts", "global_pvalue",
                       "empirical_fdr", "thresholds"))
  expect_length(null$null_counts, 15)
  # real codiversification: observed exceeds the permuted null
  expect_gt(null$observed, mean(null$null_counts))
  expect_lte(null$global_pvalue, 0.1)
  expect_gte(null$empirical_fdr, 0)
  expect_lte(null$empirical_fdr, 1)
})

test_that("codiv_null_scans preserves the global RNG", {
  skip_if_not_installed("adephylo")
  set.seed(2)
  host <- ape::rcoal(10); host$tip.label <- paste0("H", 1:10)
  sym <- host; sym$tip.label <- paste0("S", 1:10)
  df <- data.frame(Host = paste0("H", 1:10), Symbiont = paste0("S", 1:10),
                   stringsAsFactors = FALSE)

  set.seed(99); before <- runif(3)
  set.seed(99)
  invisible(codiv_null_scans(host, sym, df, n_permutations = 3,
                             min_hosts = 4, min_symbiont_tips = 7,
                             span_fraction = 1, permutations = 99,
                             methods = "hommola", cores = 1, seed = 5,
                             verbose = FALSE))
  after <- runif(3)
  expect_equal(before, after)
})

test_that("codiv_null_scans errors on a missing statistic column", {
  skip_if_not_installed("adephylo")
  host <- ape::rcoal(10); host$tip.label <- paste0("H", 1:10)
  sym <- host; sym$tip.label <- paste0("S", 1:10)
  df <- data.frame(Host = paste0("H", 1:10), Symbiont = paste0("S", 1:10),
                   stringsAsFactors = FALSE)
  expect_error(
    codiv_null_scans(host, sym, df, n_permutations = 1,
                     statistic = "Nonexistent_col",
                     min_hosts = 4, min_symbiont_tips = 7, span_fraction = 1,
                     permutations = 99, methods = "hommola", cores = 1,
                     verbose = FALSE),
    "missing"
  )
})
