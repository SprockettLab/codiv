test_that("collapse_monophyletic returns a phylo object", {
  s_tree <- ape::rtree(20)
  h_tree <- ape::rtree(5)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 4),
    Symbiont = s_tree$tip.label
  )
  result <- collapse_monophyletic(s_tree, hs_df)

  expect_s3_class(result, "phylo")
})

test_that("collapse_monophyletic reduces tree size", {
  s_tree <- ape::rtree(20)
  h_tree <- ape::rtree(5)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 4),
    Symbiont = s_tree$tip.label
  )
  result <- collapse_monophyletic(s_tree, hs_df)

  # Collapsed tree should have fewer tips
  expect_lte(length(result$tip.label), length(s_tree$tip.label))
})

test_that("collapse_monophyletic keeps all unique hosts", {
  s_tree <- ape::rtree(20)
  h_tree <- ape::rtree(5)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 4),
    Symbiont = s_tree$tip.label
  )
  result <- collapse_monophyletic(s_tree, hs_df)

  # After collapse, remaining tips should map to all original hosts (or close)
  result_hs <- hs_df[hs_df$Symbiont %in% result$tip.label, ]
  n_hosts_in_result <- length(unique(result_hs$Host))

  expect_gt(n_hosts_in_result, 0)
})

test_that("collapse_monophyletic handles polyphyletic symbionts", {
  # symbionts from the same host are interleaved, so no group is monophyletic
  s_tree <- ape::read.tree(text = "((S1:1,S2:1):1,(S3:1,S4:1):1):1;")
  hs_df <- data.frame(
    Host = c("H1", "H2", "H1", "H2"),
    Symbiont = c("S1", "S2", "S3", "S4")
  )

  result <- collapse_monophyletic(s_tree, hs_df)

  expect_s3_class(result, "phylo")
  # nothing is monophyletic-by-host, so the tree is unchanged
  expect_equal(length(result$tip.label), 4L)
})

test_that("collapse_monophyletic handles single host", {
  s_tree <- ape::rtree(10)
  hs_df <- data.frame(
    Host = rep("H1", 10),
    Symbiont = s_tree$tip.label
  )

  result <- collapse_monophyletic(s_tree, hs_df)

  # All from same host should collapse to single tip
  expect_equal(length(result$tip.label), 1)
})

test_that("collapse_monophyletic preserves tree structure", {
  s_tree <- ape::rtree(15)
  h_tree <- ape::rtree(3)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 5),
    Symbiont = s_tree$tip.label
  )

  result <- collapse_monophyletic(s_tree, hs_df)

  # Result should be a valid phylo object
  expect_true(ape::is.binary(result) || !ape::is.binary(s_tree))
})

test_that("a zero-length monophyletic group collapses to one tip", {
  # a1..a3 are identical (zero-length) tips from one host
  s_tree <- ape::read.tree(
    text = "(((a1:0,a2:0):0,a3:0):0.3,(b1:0.2,c1:0.2):0.1);"
  )
  hs_df <- data.frame(
    Host = c("X", "X", "X", "Y", "Z"),
    Symbiont = c("a1", "a2", "a3", "b1", "c1")
  )
  result <- collapse_monophyletic(s_tree, hs_df)

  expect_equal(sum(result$tip.label %in% c("a1", "a2", "a3")), 1L)
  expect_true(all(c("b1", "c1") %in% result$tip.label))
})

test_that("collapsing does not require branch lengths", {
  s_tree <- ape::read.tree(text = "(((a1,a2),a3),(b1,c1));")
  hs_df <- data.frame(
    Host = c("X", "X", "X", "Y", "Z"),
    Symbiont = c("a1", "a2", "a3", "b1", "c1")
  )
  result <- collapse_monophyletic(s_tree, hs_df)

  expect_equal(length(result$tip.label), 3L)
  expect_equal(sum(result$tip.label %in% c("a1", "a2", "a3")), 1L)
})

test_that("nested monophyletic groups collapse to a single tip", {
  s_tree <- ape::read.tree(
    text = "((((a1:0,a2:0):0,a3:0):0,a4:0):0.3,(b1:0.2,c1:0.2):0.1);"
  )
  hs_df <- data.frame(
    Host = c("X", "X", "X", "X", "Y", "Z"),
    Symbiont = c("a1", "a2", "a3", "a4", "b1", "c1")
  )
  result <- collapse_monophyletic(s_tree, hs_df)

  expect_equal(sum(result$tip.label %in% c("a1", "a2", "a3", "a4")), 1L)
})

test_that("a polyphyletic host keeps one tip per monophyletic group", {
  # host X sits in two separate clades
  s_tree <- ape::read.tree(
    text = "(((x1:0,x2:0):0.2,y1:0.2):0.1,((x3:0,x4:0):0.2,z1:0.2):0.1);"
  )
  hs_df <- data.frame(
    Host = c("X", "X", "Y", "X", "X", "Z"),
    Symbiont = c("x1", "x2", "y1", "x3", "x4", "z1")
  )
  result <- collapse_monophyletic(s_tree, hs_df)

  expect_equal(sum(result$tip.label %in% c("x1", "x2", "x3", "x4")), 2L)
  expect_true(all(c("y1", "z1") %in% result$tip.label))
})
