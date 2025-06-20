test_that("QQ Plot", {
  tmp <- qqplot(electricity)
  expect_true(ggplot2::is_ggplot(tmp))
})
