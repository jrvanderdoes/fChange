test_that("BM p-values in White Noise Tests", {
  set.seed(123)
  b <- generate_brownian_motion(250)

  res <- portmanteau_tests(b, test = 'single', lag = 10)
  expect_equal(round(res$pvalue,7), 0.262405)
  res <- portmanteau_tests(b, test = 'single', lag = 10,block_size = 1, method = 'bootstrap')
  expect_equal(round(res$pvalue,7), 0.992)

  res <- portmanteau_tests(b, test = 'multi', lag = 10, alpha = 0.01)
  expect_equal(round(res$pvalue,7), 0.7372284)

  res <- portmanteau_tests(b, test = 'spectral', kernel = bartlett_kernel,
                           bandwidth = NULL, alpha = 0.05)
  expect_equal(round(res$pvalue,7), 0.5569972)

  res <- portmanteau_tests(b, test = 'spectral', alpha = 0.1, kernel = parzen_kernel,
                           bandwidth = adaptive_bandwidth(b, kernel=parzen_kernel))
  expect_equal(round(res$pvalue,7), 0.683296)

  res <- portmanteau_tests(b, test = 'independence', components = 3, lag = 3)
  expect_equal(round(res$pvalue,7), 0.5769329)

})

test_that("Electricity statistics in White Noise Tests", {
  set.seed(123)
  b <- electricity

  res <- portmanteau_tests(b, test = 'single', lag = 10)
  expect_equal(round(res$statistic,7), 14539596)
  res <- portmanteau_tests(b, test = 'single', lag = 1, method = 'bootstrap',block_size = 1)
  expect_equal(round(res$pvalue,7), 0)

  res <- portmanteau_tests(b, test = 'multi', lag = 10, alpha = 0.01)
  expect_equal(round(res$statistic,7), 86336179)

  res <- portmanteau_tests(b, test = 'spectral', kernel = bartlett_kernel,
                           bandwidth = NULL, alpha = 0.05)
  expect_equal(round(res$statistic,7), 16.6003364)

  res <- portmanteau_tests(b, test = 'spectral', alpha = 0.1, kernel = parzen_kernel,
                           bandwidth = adaptive_bandwidth(b, kernel=parzen_kernel))
  expect_equal(round(res$statistic,7), 34.45559)

  res <- portmanteau_tests(b, test = 'independence', components = 3, lag = 3)
  expect_equal(round(res$statistic,7), 753.08642)
})
