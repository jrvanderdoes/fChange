test_that("Kernels", {
  expect_equal(truncated_kernel(-3:3/2), c(0,1,1,1,1,1,0))
  expect_equal(bartlett_kernel(-3:3/2), c(0,0,0.5,1,0.5,0,0))
  expect_equal(parzen_kernel(-3:3/2), c(0,0,0.25,1,0.25,0,0))
  expect_equal(tukey_hanning_kernel(-3:3/2), c(0,0,0.5,1,0.5,0,0))
  expect_equal(round(quadratic_spectral_kernel(-3:3/2),4),
               c(-0.0857,0.1379,0.6869,1,0.6869,0.1379,-0.0857))
  expect_equal(round(daniell_kernel(-3:3/2),4), c(0,0,0.6366,1,0.6366,0,0))
  expect_equal(flat_top_kernel(-3:3/2), c(0,0.1,0.6,1,0.6,0.1,0))
})

test_that("Adaptive Bandwidth", {
  expect_equal(round(adaptive_bandwidth(electricity,bartlett_kernel),4), 69.3986)
})

