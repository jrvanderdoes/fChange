test_that("Changes found", {
  set.seed(123)
  res <- .change_mean(electricity, M = 250)
  expect_equal(res$pvalue, 0)
  expect_equal(res$location, 124)
})

test_that("Simple changes found", {
  set.seed(123)
  data <- data.frame(rep(1,10),rep(1,10),rep(1,10),
                     rep(10,10),rep(10,10),rep(10,10))
  res <- .change_mean(data, M = 250)
  expect_equal(res$location, 3)
})

test_that("Data type error caught", {
  expect_error(.change_mean(1:10, M = 250),
               "Data type")
})

test_that("Statistic error caught", {
  expect_error(.change_mean(electricity,statistic = 'Qn', M = 250),
               "statistic")
})
