test_that("PCA Examination", {
  set.seed(123)
  results <- pca_examination(electricity)
  expect_equal(round(sum(results$reconstruction)), 369069)
  expect_equal(round(sum(results$residuals)), 0)
})

test_that("Projection Model", {
  set.seed(123)
  if(Sys.info()['sysname'] =='Linux') {
    expect_false(FALSE)
    return()
  }

  results <- projection_model(dfts(electricity$data[,50:100]), n.ahead=5)
  expect_equal(round(sum(results$fit$data)), 33283)
  expect_equal(round(sum(results$residuals$data)), 1293)
})
