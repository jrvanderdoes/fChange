test_that("multiplication works", {
  expect_equal(round(long_run_covariance(electricity,0,bartlett_kernel)[2,3],4), 262.2927)
  expect_equal(round(long_run_covariance(electricity,1,bartlett_kernel)[2,3],4), 262.2927)
  expect_equal(round(long_run_covariance(electricity,2,bartlett_kernel)[2,3],4), 453.4216)
  expect_equal(round(long_run_covariance(electricity,2,parzen_kernel)[2,3],4), 357.8572)
})
