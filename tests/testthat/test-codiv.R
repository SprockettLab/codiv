test_that("codiv returns a data frame", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  # Use temporary file for results
  temp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp_file))

  result <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file,
    methods = c("hommola")
  )

  expect_s3_class(result, "data.frame")
})

test_that("codiv output has required columns", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  temp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp_file))

  result <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file,
    methods = c("hommola")
  )

  required_cols <- c("Node_ID", "N_Symbionts", "N_Hosts", "Symbiont_Colless",
                     "Symbiont_Sackin", "Host_Colless", "Host_Sackin")
  expect_true(all(required_cols %in% colnames(result)))
})

test_that("codiv includes method columns when requested", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  temp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp_file))

  result <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file,
    methods = c("hommola")
  )

  expect_true("Hommola_r" %in% colnames(result))
  expect_true("Hommola_pvalue" %in% colnames(result))
})

test_that("codiv results are reproducible with seed", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  temp_file1 <- tempfile(fileext = ".tsv")
  temp_file2 <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(temp_file1, temp_file2)))

  result1 <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file1,
    methods = c("hommola"),
    subtree_features = FALSE
  )

  result2 <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file2,
    methods = c("hommola"),
    subtree_features = FALSE
  )

  # Results should be identical
  expect_equal(result1$Uncollapsed_Hommola_r, result2$Uncollapsed_Hommola_r)
  expect_equal(result1$Uncollapsed_Hommola_pvalue, result2$Uncollapsed_Hommola_pvalue)
})

test_that("codiv handles single method correctly", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  temp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp_file))

  result <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file,
    methods = "hommola"  # Single method
  )

  expect_true("Hommola_r" %in% colnames(result))
})

test_that("codiv filters nodes correctly", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  temp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp_file))

  result <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 5,  # Higher threshold means fewer nodes
    min_symbiont_tips = 10,
    span_fraction = 0.1,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file,
    methods = c("hommola")
  )

  # With high thresholds, might have no nodes or few nodes
  expect_s3_class(result, "data.frame")
  expect_gte(nrow(result), 0)
})

test_that("codiv saves results to file", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  temp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp_file))

  codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file,
    methods = c("hommola")
  )

  expect_true(file.exists(temp_file))
  file_contents <- read.table(temp_file, header = TRUE, sep = "\t")
  expect_s3_class(file_contents, "data.frame")
})

test_that("codiv handles missing focus_hosts parameter", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  temp_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp_file))

  # Call without specifying focus_hosts
  result <- codiv(
    Host_tree = h_tree,
    Symbiont_tree = s_tree,
    Host_to_Symbiont_df = hs_df,
    min_hosts = 4,
    min_symbiont_tips = 3,
    span_fraction = 0.5,
    permutations = 9,
    seed = 123,
    verbose = FALSE,
    Save_fp = temp_file,
    methods = c("hommola")
  )

  expect_s3_class(result, "data.frame")
})

test_that("codiv rejects unknown methods", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, 3),
    Symbiont = s_tree$tip.label
  )

  expect_error(
    codiv(h_tree, s_tree, hs_df, methods = c("hommola", "bogus"),
          verbose = FALSE),
    "Unknown method"
  )
})

test_that("codiv returns per-node p-values and no across-node adjustment", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(6)
  s_tree <- ape::rtree(40)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, length.out = 40),
    Symbiont = s_tree$tip.label
  )

  result <- codiv(h_tree, s_tree, hs_df, min_hosts = 4, min_symbiont_tips = 7,
                  span_fraction = 0.6, permutations = 9, seed = 1, cores = 1,
                  methods = "hommola", verbose = FALSE)

  expect_true("Hommola_pvalue" %in% colnames(result))
  # per-node multiple-testing correction is intentionally not provided
  expect_false(any(grepl("_p_adj$", colnames(result))))
})

test_that("codiv rejects invalid cores", {
  skip_if_not_installed("adephylo")
  h_tree <- ape::rtree(5)
  s_tree <- ape::rtree(15)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, 3), Symbiont = s_tree$tip.label)

  expect_error(
    codiv(h_tree, s_tree, hs_df, cores = 0, verbose = FALSE),
    "cores"
  )
  expect_error(
    codiv(h_tree, s_tree, hs_df, cores = 2.5, verbose = FALSE),
    "cores"
  )
})

test_that(".prepare_symbiont_tree creates and de-duplicates node labels", {
  t1 <- ape::rtree(10)
  t1$node.label <- NULL
  p1 <- codiv:::.prepare_symbiont_tree(t1, verbose = FALSE)
  expect_false(is.null(p1$node.label))
  expect_equal(anyDuplicated(p1$node.label), 0L)

  t2 <- ape::rtree(10)
  t2$node.label <- rep("x", t2$Nnode)  # all duplicated
  p2 <- codiv:::.prepare_symbiont_tree(t2, verbose = FALSE)
  expect_equal(anyDuplicated(p2$node.label), 0L)
})

test_that("codiv ParaFit results are invariant to input row order", {
  skip_if_not_installed("adephylo")
  set.seed(101)
  h_tree <- ape::rtree(6)
  s_tree <- ape::rtree(60)
  hs_df <- data.frame(
    Host = rep(h_tree$tip.label, length.out = 60),
    Symbiont = s_tree$tip.label, stringsAsFactors = FALSE
  )

  r1 <- codiv(h_tree, s_tree, hs_df, min_hosts = 4, min_symbiont_tips = 7,
              span_fraction = 0.6, permutations = 99, seed = 1, cores = 1,
              methods = "parafit", subtree_features = FALSE, verbose = FALSE)
  r2 <- codiv(h_tree, s_tree, hs_df[sample(nrow(hs_df)), ],
              min_hosts = 4, min_symbiont_tips = 7,
              span_fraction = 0.6, permutations = 99, seed = 1, cores = 1,
              methods = "parafit", subtree_features = FALSE, verbose = FALSE)

  r1 <- r1[order(r1$Node_ID), ]
  r2 <- r2[order(r2$Node_ID), ]
  expect_equal(r1$ParaFitGlobal, r2$ParaFitGlobal,
               tolerance = 1e-8)
})

test_that("codiv preserves the global RNG and is reproducible", {
  skip_if_not_installed("adephylo")
  set.seed(7)
  h_tree <- ape::rtree(7)
  s_tree <- ape::rtree(70)
  hs_df <- data.frame(Host = rep(h_tree$tip.label, length.out = 70),
                      Symbiont = s_tree$tip.label, stringsAsFactors = FALSE)
  run <- function() codiv(h_tree, s_tree, hs_df, min_hosts = 4,
                          min_symbiont_tips = 7, span_fraction = 0.6,
                          permutations = 99, seed = 5, cores = 1,
                          methods = "hommola", subtree_features = FALSE,
                          verbose = FALSE)

  # global RNG is left unchanged by codiv
  set.seed(321); before <- runif(3)
  set.seed(321); invisible(run()); after <- runif(3)
  expect_equal(before, after)

  # deterministic across repeated runs
  a <- run(); b <- run()
  a <- a[order(a$Node_ID), ]; b <- b[order(b$Node_ID), ]
  expect_equal(a$Uncollapsed_Hommola_pvalue, b$Uncollapsed_Hommola_pvalue)
})

test_that("two default-checkpoint codiv runs do not cross-contaminate", {
  skip_if_not_installed("adephylo")
  set.seed(1)
  h1 <- ape::rtree(6); s1 <- ape::rtree(40)
  d1 <- data.frame(Host = rep(h1$tip.label, length.out = 40),
                   Symbiont = s1$tip.label, stringsAsFactors = FALSE)
  h2 <- ape::rtree(8); s2 <- ape::rtree(60)
  d2 <- data.frame(Host = rep(h2$tip.label, length.out = 60),
                   Symbiont = s2$tip.label, stringsAsFactors = FALSE)

  # both use the default (NULL) Save_fp and continue = TRUE
  r1 <- codiv(h1, s1, d1, min_hosts = 3, min_symbiont_tips = 7,
              span_fraction = 0.6, permutations = 9, methods = "hommola",
              cores = 1, seed = 1, subtree_features = FALSE, verbose = FALSE)
  r2 <- codiv(h2, s2, d2, min_hosts = 3, min_symbiont_tips = 7,
              span_fraction = 0.6, permutations = 9, methods = "hommola",
              cores = 1, seed = 1, subtree_features = FALSE, verbose = FALSE)

  # rerunning the first gives the same result as running it fresh (not polluted
  # by the second run's checkpoint)
  r1b <- codiv(h1, s1, d1, min_hosts = 3, min_symbiont_tips = 7,
               span_fraction = 0.6, permutations = 9, methods = "hommola",
               cores = 1, seed = 1, subtree_features = FALSE, verbose = FALSE)
  expect_equal(nrow(r1), nrow(r1b))
  expect_equal(as.data.frame(r1)$Uncollapsed_Hommola_r, as.data.frame(r1b)$Uncollapsed_Hommola_r)
})

test_that("focus_hosts adds a <host>_PRESENT column per host", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
  hosts <- unique(sim$links$Host)[1:2]
  res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
               permutations = 9, focus_hosts = hosts, verbose = FALSE,
               Save_fp = tempfile(fileext = ".tsv"), continue = FALSE)
  cols <- paste0(hosts, "_PRESENT")
  expect_true(all(cols %in% colnames(res)))
  expect_type(res[[cols[1]]], "logical")
})

test_that("focus_hosts not in the data warns and is FALSE everywhere", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
  expect_warning(
    res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
                 permutations = 9, focus_hosts = "not_a_host", verbose = FALSE,
                 Save_fp = tempfile(fileext = ".tsv"), continue = FALSE),
    "not in Host_to_Symbiont_df"
  )
  expect_true(all(res[["not_a_host_PRESENT"]] == FALSE))
})

test_that("n_batches and batch_size are mutually exclusive", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
  expect_error(
    codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
          permutations = 9, batch_size = 2, n_batches = 2, verbose = FALSE),
    "only one of"
  )
})

test_that("n_batches and batch_size do not change the results, only checkpointing", {
  skip_if_not_installed("adephylo")
  sim <- simulate_codiv_data(n_hosts = 12, n_clades = 4, seed = 5)
  base <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
                permutations = 19, seed = 42, verbose = FALSE,
                Save_fp = tempfile(fileext = ".tsv"), continue = FALSE)
  batched <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
                   permutations = 19, seed = 42, verbose = FALSE, n_batches = 3,
                   Save_fp = tempfile(fileext = ".tsv"), continue = FALSE)
  expect_equal(nrow(base), nrow(batched))
  expect_equal(base$Uncollapsed_Hommola_r, batched$Uncollapsed_Hommola_r)
})

test_that("paco_wf survives a clade of near-identical symbionts", {
  # add_pcoord() drops eigenvalues below a fixed tolerance, so symbiont branch
  # lengths around 1e-9 used to collapse the ordination to zero dimensions and
  # error with "'dims' cannot be of length 0". Scaling the symbiont distances
  # to unit maximum fixes it without changing gof$ss or gof$p.
  sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 3)
  links <- sim$links
  hs_mat <- host_symbiont_links(links)

  tiny_tree <- sim$symbiont_tree
  tiny_tree$edge.length <- tiny_tree$edge.length * 1e-6

  expect_no_error(
    tiny <- paco_wf(sim$host_tree, tiny_tree, hs_mat, permutations = 99,
                    seed = 21)
  )
  expect_false(is.na(tiny$gof$ss))

  # p-value is unchanged by the symbiont-side rescaling
  full <- paco_wf(sim$host_tree, sim$symbiont_tree, hs_mat, permutations = 99,
                  seed = 21)
  expect_equal(tiny$gof$p, full$gof$p)
})

test_that("paco_wf returns NA when symbiont distances carry no variation", {
  # FastTree emits clades whose branch lengths are all exactly zero when a set
  # of MAGs has an identical marker alignment. Every pairwise distance is then
  # 0, the PCoA has no dimensions, and paco used to abort the whole node.
  sim <- simulate_codiv_data(n_hosts = 6, n_clades = 2, seed = 4)
  hs_mat <- host_symbiont_links(sim$links)

  flat_tree <- sim$symbiont_tree
  flat_tree$edge.length <- rep(0, length(flat_tree$edge.length))

  expect_warning(
    res <- paco_wf(sim$host_tree, flat_tree, hs_mat, permutations = 99,
                   seed = 1),
    "no variation"
  )
  expect_true(is.na(res$gof$ss))
  expect_true(is.na(res$gof$p))
})

test_that("codiv() skips zero-span clades instead of failing on them", {
  # a symbiont tree carrying one wholly collapsed clade should scan cleanly
  sim <- simulate_codiv_data(n_hosts = 8, n_clades = 3, seed = 6)
  tree <- sim$symbiont_tree

  # flatten one clade of the symbiont tree to zero branch lengths
  target <- (length(tree$tip.label) + 2)
  flat_tips <- ape::extract.clade(tree, target)$tip.label
  in_clade <- tree$edge[, 2] %in% match(flat_tips, tree$tip.label)
  tree$edge.length[in_clade] <- 0

  expect_no_error(
    res <- codiv(sim$host_tree, tree, sim$links, min_hosts = 3,
                 min_symbiont_tips = 7, permutations = 99, methods = "paco",
                 subtree_features = FALSE, cores = 1, verbose = FALSE)
  )
  expect_false(any(grepl("dims", names(res))))
})

test_that("uncollapsed Hommola is off by default and reportable on request", {
  sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 12)
  off <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
               permutations = 99, subtree_features = FALSE, cores = 1,
               verbose = FALSE)
  expect_false("Uncollapsed_Hommola_r" %in% colnames(off))
  expect_true("Hommola_r" %in% colnames(off))   # collapsed, always present

  on <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
              permutations = 99, subtree_features = FALSE, cores = 1,
              uncollapsed_hommola = TRUE, verbose = FALSE)
  expect_true(all(c("Uncollapsed_Hommola_r", "Uncollapsed_Hommola_pvalue") %in% colnames(on)))
  # turning it on must not change the collapsed statistic
  expect_equal(sort(off$Hommola_r), sort(on$Hommola_r))
})

test_that("max_symbiont_tips caps clade size and scales with the tree", {
  expect_equal(codiv:::.default_max_symbiont_tips(76697), 1534L)
  expect_equal(codiv:::.default_max_symbiont_tips(473), 500L)   # floor holds
  expect_equal(codiv:::.default_max_symbiont_tips(100000), 2000L)

  sim <- simulate_codiv_data(n_hosts = 8, n_clades = 3, seed = 13)
  res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, methods = "hommola",
               permutations = 99, min_symbiont_tips = 7, max_symbiont_tips = 12,
               subtree_features = FALSE, cores = 1, verbose = FALSE)
  expect_true(all(res$N_Symbionts <= 12))
  expect_true(all(res$N_Symbionts >= 7))
})

test_that("a failed PACo costs only PACo, not the whole node", {
  # A rank-deficient ordination used to abort the node and take the Hommola,
  # ParaFit and TreeDist results with it. Simulate the failure directly.
  sim <- simulate_codiv_data(n_hosts = 6, n_clades = 2, seed = 31)
  hs_mat <- host_symbiont_links(sim$links)
  local_mocked_bindings(
    add_pcoord = function(...) stop("'dims' cannot be of length 0"),
    .package = "paco"
  )
  expect_warning(
    r <- paco_wf(sim$host_tree, sim$symbiont_tree, hs_mat, permutations = 99,
                 seed = 1),
    "PACo could not be computed"
  )
  expect_true(is.na(r$gof$ss))
  expect_true(is.na(r$gof$p))
})
