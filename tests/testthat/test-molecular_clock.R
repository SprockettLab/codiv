make_res <- function() {
  sim <- simulate_codiv_data(n_hosts = 14, clades = list(
    list(codiversifying = TRUE, congruence = 1, n_per_host = 3)), seed = 1)
  res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, span_fraction = 0.9,
               methods = "hommola", permutations = 199, cores = 1, seed = 1,
               subtree_features = FALSE, verbose = FALSE)
  list(res = res, host_tree = sim$host_tree)
}

test_that("molecular_clock regresses symbiont divergence on host age", {
  skip_if_not_installed("adephylo")
  d <- make_res()
  mc <- molecular_clock(d$res, d$host_tree, stat_threshold = 0.75)

  expect_named(mc, c("model", "slope", "intercept", "r_squared", "p_value",
                     "slope_ci", "per_clade"))
  expect_s3_class(mc$model, "lm")
  expect_gte(nrow(mc$per_clade), 3)
  expect_true(all(c("Node_ID", "n_hosts", "host_age", "symbiont_divergence")
                  %in% colnames(mc$per_clade)))
  # idealized co-diversification: strong positive relationship
  expect_gt(mc$slope, 0)
  expect_gt(mc$r_squared, 0.5)
})

test_that("molecular_clock accepts a custom divergence function", {
  skip_if_not_installed("adephylo")
  d <- make_res()
  mc <- molecular_clock(d$res, d$host_tree, stat_threshold = 0.75,
    symbiont_divergence = function(st)
      max(ape::node.depth.edgelength(st)[seq_along(st$tip.label)]))
  expect_gt(mc$slope, 0)
})

test_that("molecular_clock errors with too few passing clades", {
  skip_if_not_installed("adephylo")
  d <- make_res()
  expect_error(
    molecular_clock(d$res, d$host_tree, stat_threshold = 0.999,
                    p_threshold = 1e-6),
    "Fewer than 3"
  )
})

test_that("molecular_clock warns on a non-ultrametric host tree", {
  skip_if_not_installed("adephylo")
  d <- make_res()
  ht <- d$host_tree
  ht$edge.length[1] <- ht$edge.length[1] * 5   # break ultrametricity
  expect_warning(molecular_clock(d$res, ht, stat_threshold = 0.75),
                 "ultrametric")
})
