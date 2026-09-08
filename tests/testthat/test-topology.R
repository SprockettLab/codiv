test_that(".host_representatives keeps the largest clade and excludes scattered hosts", {
  # H1 monophyletic pair; H2 two scattered singletons; H3/H4/H5 clean singletons
  t <- ape::read.tree(text = "(((S1,S2),(H2a,S3)),((S4,H2b),S5));")
  host_of <- c(S1 = "H1", S2 = "H1", S3 = "H3", S4 = "H4", S5 = "H5",
               H2a = "H2", H2b = "H2")
  reps <- codiv:::.host_representatives(t, host_of[t$tip.label])

  expect_equal(reps$hosts_used, 4)          # H2 excluded
  expect_false("H2" %in% host_of[reps$keep])
  expect_true("S1" %in% reps$keep || "S2" %in% reps$keep)  # H1 clade kept
  expect_equal(reps$fraction_kept, 5 / 7)   # 2 of 7 tips discarded
})

test_that("largest-clade wins with the branch-length tiebreak", {
  # H1 has two clades of equal size (2); the tighter one (shorter branches) wins
  t <- ape::read.tree(
    text = "((S1:1,S2:1):1,((S3:0.1,S4:0.1):0.1,(S5:1,S6:1):1):1);")
  host_of <- c(S1 = "H1", S2 = "H1", S3 = "H1", S4 = "H1",
               S5 = "H1", S6 = "H1")
  # all one host, two size-2 clades {S3,S4} (tight) and {S5,S6} (loose) plus
  # {S1,S2}; the largest single-host clade here is the whole thing -> monophyletic
  reps <- codiv:::.host_representatives(t, host_of[t$tip.label])
  expect_equal(reps$hosts_used, 1)          # one host, one representative
})

test_that("topology method integrates into codiv and reports diagnostics", {
  skip_if_not_installed("adephylo")
  skip_if_not_installed("TreeDist")
  set.seed(7)
  host <- ape::rcoal(14); host$tip.label <- paste0("H", 1:14)
  sym <- host; sym$tip.label <- paste0("S", 1:14)
  sym$edge.length <- sym$edge.length * runif(length(sym$edge.length), 0.7, 1.3)
  df <- data.frame(Host = paste0("H", 1:14), Symbiont = paste0("S", 1:14),
                   stringsAsFactors = FALSE)

  res <- codiv(host, sym, df, min_hosts = 4, min_symbiont_tips = 7,
               span_fraction = 0.95, permutations = 99, methods = "topology",
               cores = 1, seed = 42, subtree_features = FALSE, verbose = FALSE)

  expect_true(all(c("Topology_congruence", "Topology_pvalue",
                    "Topology_hosts_used", "Topology_fraction_kept")
                  %in% colnames(res)))
  # identical topology (branch-length noise only) -> perfect congruence
  expect_true(all(abs(res$Topology_congruence - 1) < 1e-6))
  expect_true(all(res$Topology_fraction_kept == 1))
  # summary recognizes the method
  sm <- summary(res)
  expect_true("topology" %in% sm$method)
})

test_that("topology is opt-in, not a default method", {
  expect_false("topology" %in% eval(formals(codiv)$methods))
})
