test_that("simulate_codiv_data builds the requested clades with truth", {
  sim <- simulate_codiv_data(
    n_hosts = 10,
    clades = list(
      list(codiversifying = TRUE,  congruence = 1.0, n_per_host = 1),
      list(codiversifying = TRUE,  congruence = 0.6, n_per_host = 4),
      list(codiversifying = FALSE, n_per_host = 1),
      list(codiversifying = FALSE, n_per_host = 4)
    ), seed = 1)

  expect_named(sim, c("host_tree", "symbiont_tree", "links", "truth"))
  expect_s3_class(sim$host_tree, "phylo")
  expect_s3_class(sim$symbiont_tree, "phylo")
  expect_equal(length(sim$host_tree$tip.label), 10)
  # tips: 10*1 + 10*4 + 10*1 + 10*4 = 100
  expect_equal(length(sim$symbiont_tree$tip.label), 100)
  expect_equal(nrow(sim$links), 100)
  # every symbiont maps to a host tip, uniquely
  expect_true(all(sim$links$Host %in% sim$host_tree$tip.label))
  expect_equal(anyDuplicated(sim$links$Symbiont), 0L)
  # ground truth: two co-diversifying clades
  clade_truth <- tapply(sim$truth$codiversifying, sim$truth$clade, `[`, 1)
  expect_equal(sum(clade_truth), 2L)
})

test_that("simulate_codiv_data is reproducible and preserves the global RNG", {
  set.seed(123); before <- runif(3)
  set.seed(123)
  a <- simulate_codiv_data(n_hosts = 8, n_clades = 3, seed = 5)
  after <- runif(3)
  expect_equal(before, after)                 # global RNG untouched
  b <- simulate_codiv_data(n_hosts = 8, n_clades = 3, seed = 5)
  expect_equal(a$links, b$links)              # reproducible
})

test_that("simulate_codiv_data random mode respects n_clades", {
  sim <- simulate_codiv_data(n_hosts = 12, n_clades = 6, seed = 7)
  expect_equal(length(unique(sim$truth$clade)), 6)
})

test_that("simulate_codiv_data errors without clades or n_clades", {
  expect_error(simulate_codiv_data(n_hosts = 8), "clades")
})

test_that("simulated data runs through codiv", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 2)
  res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, min_hosts = 3,
               min_symbiont_tips = 7, span_fraction = 0.6, permutations = 99,
               methods = "hommola", cores = 1, seed = 1,
               subtree_features = FALSE, verbose = FALSE)
  expect_s3_class(res, "codiv")
})
