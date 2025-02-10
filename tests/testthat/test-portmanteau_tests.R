test_that("BM p-values in White Noise Tests", {
  set.seed(123)
  b <- generate_brownian_motion(250)

  res <- portmanteau_tests(b, test = 'single-lag', lag = 10)
  expect_equal(round(res$p_value,7), 0.2316695)
  res <- portmanteau_tests(b, test = 'single-lag', lag = 10,block_size = 1, method = 'bootstrap')
  expect_equal(round(res$p_value,7), 0.996)

  res <- portmanteau_tests(b, test = 'multi-lag', lag = 10, alpha = 0.01)
  expect_equal(round(res$p_value,7), 0.8909353)

  res <- portmanteau_tests(b, test = 'spectral', kernel = 'Bartlett', bandwidth = 'static', alpha = 0.05)
  expect_equal(round(res$p_value,7), 0.4464877)

  res <- portmanteau_tests(b, test = 'spectral', alpha = 0.1, kernel = 'Parzen', bandwidth = 'adaptive')
  expect_equal(round(res$p_value,7), 0.5009694)

  res <- portmanteau_tests(b, test = 'independence', components = 3, lag = 3)
  expect_equal(round(res$p_value,7), 0.742398)

})

test_that("Electricity statistics in White Noise Tests", {
  set.seed(123)
  b <- electricity

  res <- portmanteau_tests(b, test = 'single-lag', lag = 10)
  expect_equal(round(res$statistic,7), 14539596)
  res <- portmanteau_tests(b, test = 'single-lag', lag = 1, method = 'bootstrap',block_size = 1)
  expect_equal(round(res$p_value,7), 0)

  res <- portmanteau_tests(b, test = 'multi-lag', lag = 10, alpha = 0.01)
  expect_equal(round(res$statistic,7), 86336179)

  res <- portmanteau_tests(b, test = 'spectral', kernel = 'Bartlett', bandwidth = 'static', alpha = 0.05)
  expect_equal(round(res$statistic,7), 16.6003364)

  res <- portmanteau_tests(b, test = 'spectral', alpha = 0.1, kernel = 'Parzen', bandwidth = 'adaptive')
  expect_equal(round(res$statistic,7), 34.45559)

  res <- portmanteau_tests(b, test = 'independence', components = 3, lag = 3)
  expect_equal(round(res$statistic,7), 753.08642)
})

test_that("Check Low discrapancy if possible", {
  if(requireNamespace('fOptions',quietly = TRUE)){
    set.seed(123)
    stat <- .single_lag_test(electricity, lag = 1, method='lowdiscrepancy')$statistic
    expect_equal(stat, 14539596)

    stat <- .multi_lag_test(data = electricity, lag = 10, method='lowdiscrepancy')$statistic
    expect_equal(stat, 86336179)
  }
})
