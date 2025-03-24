test_that("Stationarity Test", {
  set.seed(1234)
  bb <- generate_brownian_motion(100,v=seq(0,1,length.out=20))
  res <- stationarity_test(bb, critical ='resample', statistic='Mn')
  expect_equal(res$pvalue, 0.431)

  res <- stationarity_test(electricity)
  expect_equal(res$pvalue, 0)
})
