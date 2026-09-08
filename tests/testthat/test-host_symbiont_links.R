test_that("host_symbiont_links creates correct binary matrix", {
  hs_df <- data.frame(
    Host = c("H1", "H1", "H2", "H2", "H3"),
    Symbiont = c("S1", "S2", "S3", "S4", "S5")
  )
  result <- host_symbiont_links(hs_df)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)  # 3 unique hosts
  expect_equal(ncol(result), 5)  # 5 unique symbionts
  expect_equal(rownames(result), c("H1", "H2", "H3"))
  expect_equal(colnames(result), c("S1", "S2", "S3", "S4", "S5"))
})

test_that("host_symbiont_links has correct associations", {
  hs_df <- data.frame(
    Host = c("H1", "H2"),
    Symbiont = c("S1", "S2")
  )
  result <- host_symbiont_links(hs_df)

  expect_equal(result["H1", "S1"], 1)
  expect_equal(result["H1", "S2"], 0)
  expect_equal(result["H2", "S1"], 0)
  expect_equal(result["H2", "S2"], 1)
})

test_that("host_symbiont_links handles multiple symbionts per host", {
  hs_df <- data.frame(
    Host = c("H1", "H1", "H1"),
    Symbiont = c("S1", "S2", "S3")
  )
  result <- host_symbiont_links(hs_df)

  expect_equal(sum(result["H1", ]), 3)  # H1 has 3 symbionts
})

test_that("host_symbiont_links returns all zeros and ones", {
  hs_df <- data.frame(
    Host = rep(c("H1", "H2"), 5),
    Symbiont = paste0("S", 1:10)
  )
  result <- host_symbiont_links(hs_df)

  expect_true(all(result %in% c(0, 1)))
})

test_that("host_symbiont_links warns on duplicate symbionts", {
  hs_df <- data.frame(
    Host = c("H1", "H2"),
    Symbiont = c("S1", "S1")  # Duplicate symbiont
  )

  expect_warning(
    host_symbiont_links(hs_df),
    "duplicated"
  )
})

test_that("host_symbiont_links handles single association", {
  hs_df <- data.frame(Host = "H1", Symbiont = "S1")
  result <- host_symbiont_links(hs_df)

  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 1)
  expect_equal(result[1, 1], 1)
})
