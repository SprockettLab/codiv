test_that("host_codiv_summary ranks hosts in co-diversifying clades", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 10, clades = list(
    list(codiversifying = TRUE,  congruence = 1.0, n_per_host = 2),
    list(codiversifying = FALSE, n_per_host = 2)), seed = 3)
  res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, min_hosts = 3,
               min_symbiont_tips = 7, span_fraction = 0.6, permutations = 199,
               methods = "hommola", cores = 1, seed = 1,
               subtree_features = FALSE, verbose = FALSE)

  hs <- host_codiv_summary(res, stat_threshold = 0.75, p_threshold = 0.01)
  expect_named(hs, c("Host", "n_codiversifying_nodes", "n_scanned_nodes",
                     "prop_codiversifying", "max_statistic"))
  # the co-diversifying clade spans all hosts, so all should appear
  expect_setequal(hs$Host, sim$host_tree$tip.label)
  expect_true(all(hs$n_codiversifying_nodes > 0))
  expect_true(all(hs$prop_codiversifying <= 1))
  # ordered by n_codiversifying_nodes descending
  expect_false(is.unsorted(rev(hs$n_codiversifying_nodes)))
})

test_that("host_codiv_summary finds few co-diversifying hosts under the null", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 10, clades = replicate(3,
    list(codiversifying = FALSE, n_per_host = 1), simplify = FALSE), seed = 11)
  res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, min_hosts = 3,
               min_symbiont_tips = 7, span_fraction = 1, permutations = 199,
               methods = "hommola", cores = 1, seed = 1,
               subtree_features = FALSE, verbose = FALSE)
  hs <- host_codiv_summary(res, stat_threshold = 0.75, p_threshold = 0.01)
  # random data: at most a handful of co-diversifying nodes across all hosts
  expect_lt(sum(hs$n_codiversifying_nodes), sum(hs$n_scanned_nodes))
})

test_that("host_codiv_summary errors on a missing column", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 1)
  res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, min_hosts = 3,
               min_symbiont_tips = 7, span_fraction = 1, permutations = 99,
               methods = "hommola", cores = 1, seed = 1,
               subtree_features = FALSE, verbose = FALSE)
  expect_error(host_codiv_summary(res, statistic = "Nonexistent"), "missing")
})
