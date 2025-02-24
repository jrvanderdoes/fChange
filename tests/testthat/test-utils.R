test_that("Specify Decimal Check", {
  test <- .specify_decimal(1,2)
  expect_equal(test, '1.00')
  expect_equal(nchar(test), 4)

  test <- .specify_decimal(1.2345,1)
  expect_equal(test, '1.2')
  expect_equal(nchar(test), 3)
})

test_that("Select n Check", {
  expect_equal(.select_n(1:10,1), 1)
  expect_equal(.select_n(1:10,2), c(1, 10))
  expect_equal(.select_n(1:10,6), c(1, 3, 5, 6, 8, 10))
})

test_that("Get Chunks Check", {
  expect_equal(length(.getChunks(1:100, 1)), 100)
  expect_equal(length(.getChunks(1:100, 2)), 2)
  expect_equal(length(.getChunks(1:100, 4)[[1]]), 25)
})

test_that("Bootstrap Check", {
  set.seed(123)
  expect_equal(.bootstrap(electricity,1,M=2)[[2]][2,2], 29.7)
  expect_equal(length(.bootstrap(electricity,1,M=10)), 10)
  expect_equal(length(.bootstrap(electricity,1,M=10,type = 'separate')), 10)
  expect_equal(
    round(.bootstrap(electricity,1,M=10,type = 'separate',
                     fn = .characteristic_statistic,statistic='Tn')[[2]],4), 0.1249)
})
