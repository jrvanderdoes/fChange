test_that("PCA Examination", {
  results <- pca_examination(electricity)
  expect_equal(round(sum(results$reconstruction)), 369069)
  expect_equal(round(sum(results$residuals)), 0)
})

test_that("Projection Model", {
  result <- projection_model(dfts(electricity$data[,50:150]), n.ahead=10)
  expect_equal(round(sum(results$reconstruction)), 369069)
  expect_equal(round(sum(results$residuals)), 0)
})
