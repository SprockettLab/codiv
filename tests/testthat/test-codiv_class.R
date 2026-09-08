make_result <- function(uncollapsed_hommola = FALSE) {
  skip_if_not_installed("adephylo")
  set.seed(42)
  h <- ape::rtree(8); s <- ape::rtree(120)
  hs <- data.frame(Host = rep(h$tip.label, 15), Symbiont = s$tip.label,
                   stringsAsFactors = FALSE)
  codiv(h, s, hs, min_hosts = 4, min_symbiont_tips = 7, span_fraction = 0.5,
        permutations = 99, seed = 123, cores = 1,
        methods = c("hommola", "paco"),
        uncollapsed_hommola = uncollapsed_hommola, verbose = FALSE)
}

test_that("codiv returns a codiv object that is also a data frame", {
  res <- make_result()
  expect_s3_class(res, "codiv")
  expect_s3_class(res, "data.frame")
  expect_false(is.null(attr(res, "codiv_params")))
  expect_equal(attr(res, "codiv_params")$permutations, 99)
})

test_that("print.codiv returns its input invisibly", {
  res <- make_result()
  expect_output(print(res), "codiversification scan")
  expect_invisible(print(res))
})

test_that("summary.codiv returns a row per method plus collapsed Hommola", {
  # the second Hommola row only exists when the uncollapsed test was run
  res <- make_result(uncollapsed_hommola = TRUE)
  sm <- summary(res)
  expect_s3_class(sm, "data.frame")
  expect_setequal(sm$method, c("hommola", "hommola_uncollapsed", "paco"))
  # the collapsed test is the primary row; the opt-in one is supplementary
  expect_equal(sm$statistic[sm$method == "hommola"], "Hommola_r")
  expect_equal(sm$statistic[sm$method == "hommola_uncollapsed"],
               "Uncollapsed_Hommola_r")
  expect_true(all(c("n_sig_p", "statistic") %in% colnames(sm)))
})

test_that("summary.codiv reports only the collapsed Hommola by default", {
  sm <- summary(make_result())
  expect_setequal(sm$method, c("hommola", "paco"))
  expect_equal(sm$statistic[sm$method == "hommola"], "Hommola_r")
})

test_that("as.data.frame.codiv strips the class and attribute", {
  res <- make_result()
  d <- as.data.frame(res)
  expect_identical(class(d), "data.frame")
  expect_null(attr(d, "codiv_params"))
})

test_that("codiv object still supports data-frame subsetting", {
  res <- make_result()
  expect_true(nrow(res) >= 1)
  sub <- res[res$N_Hosts >= 4, ]
  expect_true(is.data.frame(sub))
  expect_true("Hommola_r" %in% colnames(res))
})
