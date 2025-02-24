test_that("QQ Plot", {
  tmp <- qqplot(electricity)
  expect_equal(class(tmp)[1], 'gg')
})
