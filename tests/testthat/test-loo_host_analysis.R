test_that("loo_host_analysis returns the expected structure", {
  skip_if_not_installed("adephylo")
  set.seed(42)
  h <- ape::rtree(6); s <- ape::rtree(120)
  hs <- data.frame(Host = rep(h$tip.label, length.out = 120),
                   Symbiont = s$tip.label, stringsAsFactors = FALSE)

  loo <- loo_host_analysis(h, s, hs, min_hosts = 4, min_symbiont_tips = 7,
                           span_fraction = 0.6, permutations = 99,
                           methods = c("hommola", "paco"), cores = 1,
                           seed = 123, verbose = FALSE)

  expect_named(loo, c("full_results", "loo_results", "host_importance"))
  expect_s3_class(loo$full_results, "codiv")
  expect_length(loo$loo_results, 6)
  expect_true(all(vapply(loo$loo_results,
                         function(r) inherits(r, "codiv"), logical(1))))
  # one row per host, ordered by impact_score (descending)
  expect_equal(nrow(loo$host_importance), 6)
  expect_setequal(loo$host_importance$host, h$tip.label)
  expect_true(!is.unsorted(rev(loo$host_importance$impact_score)))
  expect_true(all(c("hommola_delta_n_sig", "paco_delta_n_sig", "impact_score")
                  %in% colnames(loo$host_importance)))
})

test_that("loo runs genuinely recompute with the host removed", {
  skip_if_not_installed("adephylo")
  set.seed(42)
  h <- ape::rtree(6); s <- ape::rtree(120)
  hs <- data.frame(Host = rep(h$tip.label, length.out = 120),
                   Symbiont = s$tip.label, stringsAsFactors = FALSE)

  loo <- loo_host_analysis(h, s, hs, hosts_to_test = c("t1", "t2"),
                           min_hosts = 4, min_symbiont_tips = 7,
                           span_fraction = 0.6, permutations = 99,
                           methods = "hommola", cores = 1, seed = 123,
                           verbose = FALSE)

  expect_length(loo$loo_results, 2)
  # a host-removed run produces valid (non all-NA) statistics
  expect_false(all(is.na(loo$loo_results[["t1"]]$Hommola_r)))
  # fewer or equal nodes than the full run (removing a host + its symbionts)
  expect_lte(nrow(loo$loo_results[["t1"]]), nrow(loo$full_results))
})

test_that("loo_host_analysis rejects unknown hosts", {
  skip_if_not_installed("adephylo")
  h <- ape::rtree(6); s <- ape::rtree(30)
  hs <- data.frame(Host = rep(h$tip.label, 5), Symbiont = s$tip.label,
                   stringsAsFactors = FALSE)
  expect_error(
    loo_host_analysis(h, s, hs, hosts_to_test = "not_a_host", verbose = FALSE),
    "not present"
  )
})

test_that("codiv drops symbionts with no host association", {
  skip_if_not_installed("adephylo")
  set.seed(1)
  h <- ape::rtree(5); s <- ape::rtree(40)
  # only associate the first 30 symbionts; 10 tips are orphans
  hs <- data.frame(Host = rep(h$tip.label, length.out = 30),
                   Symbiont = s$tip.label[1:30], stringsAsFactors = FALSE)
  expect_message(
    res <- codiv(h, s, hs, min_hosts = 4, min_symbiont_tips = 7,
                 span_fraction = 0.6, permutations = 99, methods = "hommola",
                 cores = 1, seed = 1, subtree_features = FALSE),
    "no host association"
  )
  expect_false(all(is.na(res$Hommola_r)))
})
