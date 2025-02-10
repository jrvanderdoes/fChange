test_that("Imputation", {
  set.seed(123)
  temp <- data.frame(c(NA,NA,3:9,NA),
                     c(NA,stats::rnorm(2),NA, stats::rnorm(6)),
                     stats::rnorm(10),
                     c(stats::rnorm(4),rep(NA,3), stats::rnorm(3)),
                     rep(NA,10),
                     c(stats::rnorm(1), rep(NA,9)),
                     c(stats::rnorm(9),NA),
                     stats::rnorm(10),
                     stats::rnorm(10),
                     c(NA,NA,3:9, NA))
  temp_imp <- impute(X = temp, method='zero')
  expect_equal(temp_imp$data[2,1],0)
  expect_equal(temp_imp$data[3,5],0)


  temp_imp <- impute(X = temp, method='mean_obs')
  expect_equal(temp_imp$data[2,1],6)
  expect_equal(temp_imp$data[3,5],as.numeric(NA))

  temp_imp <- impute(X = temp, method='median_obs')
  expect_equal(temp_imp$data[2,1],6)
  expect_equal(temp_imp$data[3,5],as.numeric(NA))


  temp_imp <- impute(X = temp, method='mean_data')
  expect_equal(round(temp_imp$data[2,1],4),-0.1958)
  expect_equal(round(temp_imp$data[3,5],4),0.5324)

  temp_imp <- impute(X = temp, method='median_data')
  expect_equal(round(temp_imp$data[2,1],4),-0.4243)
  expect_equal(round(temp_imp$data[3,5],4),-0.1460)


  temp_imp <- impute(temp, method='linear', obs_share_data=TRUE)
  expect_equal(temp_imp$data[2,1],2)
  expect_equal(round(temp_imp$data[4,2],4),0.6643)
  expect_equal(round(temp_imp$data[3,5],4),-0.9146)

  expect_error(impute(temp, method='functional'))
  temp_imp <- impute(temp[,3:4], method='functional')
  expect_equal(round(temp_imp$data[5,2],4),-0.9702)
})
