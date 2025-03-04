test_that("Confidence Interval", {
  set.seed(1234)
  X <- cbind(generate_brownian_motion(200,v=seq(0,1,0.05))$data,
             generate_brownian_motion(100,v=seq(0,1,0.05))$data+0.1,
             generate_brownian_motion(150,v=seq(0,1,0.05))$data-0.05)
  res <- confidence_interval(X,c(200,300))
  expect_equal(res$change,c(200, 300))
  expect_equal(round(res$lower,4),c(44.2607, 229.7510))
  expect_equal(round(res$upper,4),c(374.7260, 378.8133))

  expect_error(confidence_interval(X,c()))
  expect_error(confidence_interval(X,c(0)))
})
