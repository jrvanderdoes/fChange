test_that("Check KPSS Test", {
  set.seed(123)
  res <- res1 <- rep(NA,2)

  dat <- generate_brownian_motion(100,v=seq(0,1,length.out=20))
  res1[1] <- kpss_test(dat)$pvalue
  res1[2] <- kpss_test(dat, method="resample")$pvalue

  for(j in 2:100){
    dat$data[,j] <- dat$data[,j]+dat$data[,j-1]
  }
  res[1] <- kpss_test(dat)$pvalue
  res[2] <- kpss_test(dat, method="resample")$pvalue

  expect_equal(res1, c(0.527, 0.504))
  expect_equal(res, c(0.000, 0.000))
})
