test_that("multiplication works", {
  tmp <- qqplot(funts(electricity))
  expect_equal(class(tmp)[1], 'gg')
})
