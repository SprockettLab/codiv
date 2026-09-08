test_that("check_inputs accepts valid inputs", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 4),
    Symbiont = s_tree$tip.label
  )

  expect_invisible(
    check_inputs(h_tree, s_tree, hs_df, min_hosts = 4, min_symbiont_tips = 7,
                 span_fraction = 0.5, permutations = 99)
  )
})

test_that("check_inputs rejects non-phylo Host_tree", {
  h_tree <- "not a tree"
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = c("h1", "h2"), Symbiont = c("s1", "s2"))

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "must be phylo objects"
  )
})

test_that("check_inputs rejects non-phylo Symbiont_tree", {
  h_tree <- ape::rtree(5)
  s_tree <- list(tip.label = c("s1", "s2"))  # Not phylo
  hs_df <- data.frame(Host = c("h1"), Symbiont = c("s1"))

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "must be phylo objects"
  )
})

test_that("check_inputs rejects min_hosts < 3", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, min_hosts = 2, min_symbiont_tips = 7,
                 span_fraction = 0.5, permutations = 99),
    ">= 3"
  )
})

test_that("check_inputs rejects span_fraction outside (0,1]", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, span_fraction = 0, permutations = 99),
    "between 0"
  )
  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, span_fraction = 1.5, permutations = 99),
    "between 0"
  )
})

test_that("check_inputs warns on low min_symbiont_tips", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)

  expect_warning(
    check_inputs(h_tree, s_tree, hs_df, 4, min_symbiont_tips = 5,
                 span_fraction = 0.5, permutations = 99),
    "7 or more"
  )
})

test_that("check_inputs warns on low permutations", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)

  expect_warning(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, span_fraction = 0.5,
                 permutations = 50),
    "99 or more"
  )
})

test_that("check_inputs rejects missing Host column", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Symbiont = s_tree$tip.label[1:5])  # Missing Host column

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "'Host' and 'Symbiont'"
  )
})

test_that("check_inputs rejects missing Symbiont column", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 2))  # Missing Symbiont column

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "'Host' and 'Symbiont'"
  )
})

test_that("check_inputs rejects symbionts not in Symbiont_tree", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(
    Host = c("h1", "h2"),
    Symbiont = c("fake_symbiont_1", "fake_symbiont_2")
  )

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "not tips in"
  )
})

test_that("check_inputs rejects hosts not in Host_tree", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(
    Host = c("fake_host_1", "fake_host_2"),
    Symbiont = s_tree$tip.label[1:2]
  )

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "not tips in"
  )
})

test_that("check_inputs rejects trees without branch lengths", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  s_tree$edge.length <- NULL  # strip branch lengths
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "branch lengths"
  )
})

test_that("check_inputs rejects duplicated tip labels", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  s_tree$tip.label[2] <- s_tree$tip.label[1]  # create a duplicate
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "duplicated tip labels"
  )
})

test_that("check_inputs rejects NA in the linkage data frame", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)
  hs_df$Host[1] <- NA

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "must not contain NA"
  )
})

test_that("check_inputs warns on non-bifurcating trees", {
  h_tree <- ape::rtree(5)
  # polytomies inside the tree but a bifurcating root, so this exercises the
  # bifurcation warning rather than the unrooted-tree error
  s_tree <- ape::read.tree(
    text = "(((a:1,b:1,c:1):1,(d:1,e:1):1):1,(f:1,g:1,h:1):1);")
  hs_df <- data.frame(Host = rep(h_tree$tip.label, length.out = 8),
                      Symbiont = s_tree$tip.label)

  expect_warning(
    check_inputs(h_tree, s_tree, hs_df, 4, 7, 0.5, 99),
    "bifurcating"
  )
})

test_that("check_inputs rejects non-integer count parameters", {
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(20)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 4), Symbiont = s_tree$tip.label)

  expect_error(
    check_inputs(h_tree, s_tree, hs_df, min_hosts = 4.5, min_symbiont_tips = 7,
                 span_fraction = 0.5, permutations = 99),
    "positive integer"
  )
  expect_error(
    check_inputs(h_tree, s_tree, hs_df, min_hosts = 4, min_symbiont_tips = 7,
                 span_fraction = "half", permutations = 99),
    "span_fraction"
  )
})

test_that("an unrooted symbiont tree is rejected", {
  set.seed(20)
  host <- ape::rtree(6)
  sym <- ape::rtree(30, rooted = FALSE)
  links <- data.frame(Host = rep(host$tip.label, 5), Symbiont = sym$tip.label)
  expect_error(check_inputs(host, sym, links, 3, 7, 0.25, 99), "unrooted")
  expect_error(check_inputs(host, sym, links, 3, 7, 0.25, 99),
               "root_at_midpoint")
})

test_that("multi2di on an unrooted tree is flagged by its zero-length root", {
  # this passes is.rooted() but nobody rooted it
  set.seed(21)
  host <- ape::rtree(6)
  sym <- ape::multi2di(ape::rtree(30, rooted = FALSE))
  links <- data.frame(Host = rep(host$tip.label, 5), Symbiont = sym$tip.label)
  expect_true(ape::is.rooted(sym))
  expect_warning(check_inputs(host, sym, links, 3, 7, 0.25, 99),
                 "zero-length branch at the root")
})

test_that("a properly rooted tree passes both rooting checks", {
  set.seed(22)
  host <- ape::rtree(6)
  sym <- ape::multi2di(castor::root_at_midpoint(ape::rtree(30, rooted = FALSE)))
  links <- data.frame(Host = rep(host$tip.label, 5), Symbiont = sym$tip.label)
  expect_silent(check_inputs(host, sym, links, 3, 7, 0.25, 99))
})

test_that("an unrooted host tree only warns", {
  set.seed(23)
  host <- ape::rtree(6, rooted = FALSE)
  sym <- ape::rtree(30)
  links <- data.frame(Host = rep(host$tip.label, 5), Symbiont = sym$tip.label)
  expect_warning(check_inputs(host, sym, links, 3, 7, 0.25, 99),
                 "Host_tree` is unrooted")
})
