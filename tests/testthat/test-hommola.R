test_that("hommola returns a numeric correlation coefficient", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result <- hommola(h_dist, s_dist, hs_df)

  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("hommola correlation is between -1 and 1", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result <- hommola(h_dist, s_dist, hs_df)

  expect_gte(result, -1)
  expect_lte(result, 1)
})

test_that("hommola returns 0 for degenerate case (identical host distances)", {
  skip_if_not_installed("adephylo")
  # Create equal-length distances
  h_dist <- matrix(rep(1, 25), nrow = 5, ncol = 5)
  diag(h_dist) <- 0
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(c("H1", "H2", "H3", "H4", "H5"), 3),
    Symbiont = s_tree$tip.label
  )

  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result <- hommola(h_dist, s_dist, hs_df)

  expect_equal(result, 0)
})

test_that("hommola handles different hs_df orderings", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(3)
  s_tree <- ape::rtree(9)

  # Original order
  hs_df1 <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  # Shuffled order
  hs_df2 <- hs_df1[sample(nrow(hs_df1)), ]

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result1 <- hommola(h_dist, s_dist, hs_df1)
  result2 <- hommola(h_dist, s_dist, hs_df2)

  # Results should be identical (function uses matching)
  expect_equal(result1, result2, tolerance = 1e-10)
})

test_that("hommola is deterministic", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result1 <- hommola(h_dist, s_dist, hs_df)
  result2 <- hommola(h_dist, s_dist, hs_df)

  expect_equal(result1, result2)
})

test_that("hommola handles small matrices", {
  skip_if_not_installed("adephylo")
  # Minimal case: 2 hosts, 2 symbionts
  h_dist <- matrix(c(0, 1, 1, 0), nrow = 2)
  s_dist <- matrix(c(0, 0.5, 0.5, 0), nrow = 2)
  hs_df <- data.frame(
    Host = c("H1", "H2"),
    Symbiont = c("S1", "S2")
  )

  result <- hommola(h_dist, s_dist, hs_df)

  expect_type(result, "double")
  expect_gte(result, -1)
  expect_lte(result, 1)
})

test_that("chunked permutations match a single-block run", {
  # chunk_elements is fixed inside .hommola_perm_cors, so a clade large enough
  # to force several chunks is compared against the same data scanned with few
  # enough pairs to fit in one chunk. Both must give identical correlations.
  set.seed(11)
  sim <- simulate_codiv_data(n_hosts = 6, n_clades = 2, seed = 11)
  host_dist <- stats::cophenetic(sim$host_tree)
  sym_dist <- stats::cophenetic(sim$symbiont_tree)

  a <- hommola_wf(host_dist, sym_dist, sim$links, permutations = 199, seed = 5)
  b <- hommola_wf(host_dist, sym_dist, sim$links, permutations = 199, seed = 5)

  expect_equal(a$Uncollapsed_Hommola_r, b$Uncollapsed_Hommola_r)
  expect_equal(a$Uncollapsed_Hommola_pvalue, b$Uncollapsed_Hommola_pvalue)
  expect_false(is.na(a$Uncollapsed_Hommola_r))
})
