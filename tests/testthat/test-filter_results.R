# small helper: build a fake codiv result from an explicit nesting structure
fake_results <- function() {
  # ((a,b),(c,d)) plus a disjoint (e,f); calls are clades of that tree
  tips <- list(
    big    = c("a", "b", "c", "d"),
    left   = c("a", "b"),
    right  = c("c", "d"),
    other  = c("e", "f")
  )
  data.frame(
    Node_ID = names(tips),
    Symbiont_Tree = vapply(tips, function(x)
      paste0("(", paste(x, collapse = ","), ");"), character(1)),
    Uncollapsed_Hommola_r = c(0.5, 0.9, 0.4, 0.8),
    Uncollapsed_Hommola_pvalue = c(0.001, 0.001, 0.20, 0.001),
    Hommola_r = c(0.5, 0.9, 0.4, 0.8),
    Homo_sapiens_PRESENT = c(TRUE, TRUE, FALSE, TRUE),
    Chlorocebus_sabaeus_PRESENT = c(TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("column expressions subset like subset()", {
  r <- filter_results(fake_results(), Uncollapsed_Hommola_pvalue < 0.05)
  expect_equal(sort(r$Node_ID), c("big", "left", "other"))
})

test_that("multiple expressions combine with AND", {
  r <- filter_results(fake_results(), Uncollapsed_Hommola_pvalue < 0.05, Uncollapsed_Hommola_r > 0.85)
  expect_equal(r$Node_ID, "left")
})

test_that("NA in a filtering expression drops the row", {
  df <- fake_results()
  df$Uncollapsed_Hommola_r[1] <- NA
  r <- filter_results(df, Uncollapsed_Hommola_r > 0.1)
  expect_false("big" %in% r$Node_ID)
})

test_that("hosts filters on the _PRESENT columns", {
  df <- fake_results()
  expect_equal(sort(filter_results(df, hosts = "Homo_sapiens")$Node_ID),
               c("big", "left", "other"))
  expect_equal(filter_results(df, hosts = c("Homo_sapiens", "Chlorocebus_sabaeus"))$Node_ID,
               c("big", "other"))
  expect_equal(sort(filter_results(df, hosts = c("Homo_sapiens", "Chlorocebus_sabaeus"),
                                   hosts_match = "any")$Node_ID),
               c("big", "left", "other", "right"))
})

test_that("an untracked host gives an informative error", {
  expect_error(filter_results(fake_results(), hosts = "Pan_troglodytes"),
               "focus_hosts")
  expect_error(filter_results(fake_results(), hosts = "Pan_troglodytes"),
               "Homo_sapiens")
})

test_that("nested = keep_all is the default and changes nothing", {
  df <- fake_results()
  expect_equal(nrow(filter_results(df)), 4)
  expect_equal(filter_results(df)$Node_ID, filter_results(df, nested = "keep_all")$Node_ID)
})

test_that("outermost keeps the largest call per nested series", {
  r <- filter_results(fake_results(), nested = "outermost")
  expect_equal(sort(r$Node_ID), c("big", "other"))
})

test_that("innermost keeps calls that contain no other surviving call", {
  r <- filter_results(fake_results(), nested = "innermost")
  # 'big' contains left and right, so it is excluded; the rest are minimal
  expect_equal(sort(r$Node_ID), c("left", "other", "right"))
  # and here innermost returns more rows than outermost
  expect_gt(nrow(r), nrow(filter_results(fake_results(), nested = "outermost")))
})

test_that("best picks the strongest statistic in each series", {
  r <- filter_results(fake_results(), nested = "best")
  # 'left' (0.9) beats 'big' (0.5) and 'right' (0.4) in that family
  expect_equal(sort(r$Node_ID), c("left", "other"))
})

test_that("best falls back to the outermost call when the statistic is all NA", {
  df <- fake_results()
  df$Hommola_r <- NA_real_
  r <- filter_results(df, nested = "best")
  expect_equal(sort(r$Node_ID), c("big", "other"))
})

test_that("filtering runs before de-nesting", {
  # dropping 'big' on p-value must promote 'left' and 'right' to outermost,
  # rather than selecting 'big' first and then discarding the whole family
  df <- fake_results()
  df$Uncollapsed_Hommola_pvalue <- c(0.9, 0.001, 0.001, 0.001)
  r <- filter_results(df, Uncollapsed_Hommola_pvalue < 0.05, nested = "outermost")
  expect_equal(sort(r$Node_ID), c("left", "other", "right"))
})

test_that("result keeps the codiv class and records the filter", {
  r <- filter_results(fake_results(), Uncollapsed_Hommola_pvalue < 0.05, nested = "outermost")
  expect_s3_class(r, "codiv")
  f <- attr(r, "codiv_filter")
  expect_equal(f$nested, "outermost")
  expect_match(f$conditions, "Uncollapsed_Hommola_pvalue")
})

test_that("an empty result is handled without error", {
  expect_equal(nrow(filter_results(fake_results(), Uncollapsed_Hommola_r > 10,
                                   nested = "outermost")), 0)
})

test_that("a non-logical filtering expression is rejected", {
  expect_error(filter_results(fake_results(), Uncollapsed_Hommola_r + 1), "logical")
})

test_that("sort_by_filters orders by each condition's implied direction", {
  df <- fake_results()
  r <- filter_results(df, Uncollapsed_Hommola_pvalue < 0.5, Uncollapsed_Hommola_r > 0.1,
                      sort_by_filters = TRUE)
  # p ascending first, then r descending as the tiebreak
  expect_equal(r$Uncollapsed_Hommola_pvalue, sort(r$Uncollapsed_Hommola_pvalue))
  ties <- r$Uncollapsed_Hommola_pvalue == r$Uncollapsed_Hommola_pvalue[1]
  expect_equal(r$Uncollapsed_Hommola_r[ties], sort(r$Uncollapsed_Hommola_r[ties], decreasing = TRUE))
})

test_that("sort direction follows the operator", {
  df <- fake_results()
  up <- filter_results(df, Uncollapsed_Hommola_r > 0.1, sort_by_filters = TRUE)
  expect_equal(up$Uncollapsed_Hommola_r, sort(up$Uncollapsed_Hommola_r, decreasing = TRUE))
  dn <- filter_results(df, Uncollapsed_Hommola_r < 10, sort_by_filters = TRUE)
  expect_equal(dn$Uncollapsed_Hommola_r, sort(dn$Uncollapsed_Hommola_r))
})

test_that("a backwards comparison flips the sort direction", {
  df <- fake_results()
  a <- filter_results(df, Uncollapsed_Hommola_r > 0.1, sort_by_filters = TRUE)
  b <- filter_results(df, 0.1 < Uncollapsed_Hommola_r, sort_by_filters = TRUE)
  expect_equal(a$Node_ID, b$Node_ID)
})

test_that("sorting happens after de-nesting", {
  df <- fake_results()
  r <- filter_results(df, Uncollapsed_Hommola_r > 0.1, nested = "outermost",
                      sort_by_filters = TRUE)
  expect_equal(sort(r$Node_ID), c("big", "other"))
  expect_equal(r$Uncollapsed_Hommola_r, sort(r$Uncollapsed_Hommola_r, decreasing = TRUE))
})

test_that("sort_by_filters warns when no condition gives a direction", {
  expect_warning(filter_results(fake_results(), hosts = "Homo_sapiens",
                                sort_by_filters = TRUE),
                 "simple `column op constant`")
})

test_that("the applied sort is recorded on the result", {
  r <- filter_results(fake_results(), Uncollapsed_Hommola_pvalue < 0.5, Uncollapsed_Hommola_r > 0.1,
                      sort_by_filters = TRUE)
  expect_equal(attr(r, "codiv_filter")$sorted_by,
               c("Uncollapsed_Hommola_pvalue (asc)", "Uncollapsed_Hommola_r (desc)"))
})
