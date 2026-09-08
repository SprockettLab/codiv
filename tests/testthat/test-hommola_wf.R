test_that("hommola_wf returns list with Uncollapsed_Hommola_r and Uncollapsed_Hommola_pvalue", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result <- hommola_wf(h_dist, s_dist, hs_df, permutations = 49, seed = 123)

  expect_type(result, "list")
  expect_named(result, c("Uncollapsed_Hommola_r", "Uncollapsed_Hommola_pvalue"))
})

test_that("hommola_wf correlation is valid", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result <- hommola_wf(h_dist, s_dist, hs_df, permutations = 49, seed = 123)

  expect_gte(result$Uncollapsed_Hommola_r, -1)
  expect_lte(result$Uncollapsed_Hommola_r, 1)
})

test_that("hommola_wf p-value is between 0 and 1", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result <- hommola_wf(h_dist, s_dist, hs_df, permutations = 49, seed = 123)

  expect_gte(result$Uncollapsed_Hommola_pvalue, 0)
  expect_lte(result$Uncollapsed_Hommola_pvalue, 1)
})

test_that("hommola_wf is reproducible with seed", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result1 <- hommola_wf(h_dist, s_dist, hs_df, permutations = 49, seed = 123)
  result2 <- hommola_wf(h_dist, s_dist, hs_df, permutations = 49, seed = 123)

  expect_equal(result1$Uncollapsed_Hommola_r, result2$Uncollapsed_Hommola_r)
  expect_equal(result1$Uncollapsed_Hommola_pvalue, result2$Uncollapsed_Hommola_pvalue)
})

test_that("hommola_wf differs with different seed", {
  skip_if_not_installed("adephylo")
  # fix tree generation so the test is deterministic, and use a larger data set
  # with more permutations so two seeds reliably yield different null draws
  set.seed(20240101)
  h_tree <- ape::rtree(8)
  s_tree <- ape::rtree(48)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 6),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result1 <- hommola_wf(h_dist, s_dist, hs_df, permutations = 199, seed = 123)
  result2 <- hommola_wf(h_dist, s_dist, hs_df, permutations = 199, seed = 456)

  # observed correlation is seed-independent
  expect_equal(result1$Uncollapsed_Hommola_r, result2$Uncollapsed_Hommola_r)
  # but the permutation-based p-values should differ
  expect_false(isTRUE(all.equal(result1$Uncollapsed_Hommola_pvalue, result2$Uncollapsed_Hommola_pvalue)))
})

test_that("hommola_wf handles NA seed", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result <- hommola_wf(h_dist, s_dist, hs_df, permutations = 49, seed = NA)

  expect_type(result$Uncollapsed_Hommola_r, "double")
  expect_type(result$Uncollapsed_Hommola_pvalue, "double")
})

test_that("hommola_wf with high permutations is more accurate", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  h_dist <- adephylo::distTips(h_tree, method = "patristic")
  s_dist <- adephylo::distTips(s_tree, method = "patristic")

  result_low <- hommola_wf(h_dist, s_dist, hs_df, permutations = 9, seed = 123)
  result_high <- hommola_wf(h_dist, s_dist, hs_df, permutations = 999, seed = 123)

  # Both should have same correlation (observed value)
  expect_equal(result_low$Uncollapsed_Hommola_r, result_high$Uncollapsed_Hommola_r, tolerance = 1e-10)
  # But p-values should be different (more precision with more permutations)
  expect_false(result_low$Uncollapsed_Hommola_pvalue == result_high$Uncollapsed_Hommola_pvalue)
})
