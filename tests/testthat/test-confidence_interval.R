test_that("Confidence Interval", {
  set.seed(1234)
  X <- cbind(generate_brownian_motion(200,v=seq(0,1,0.05))$data,
             generate_brownian_motion(100,v=seq(0,1,0.05))$data+0.1,
             generate_brownian_motion(150,v=seq(0,1,0.05))$data-0.05)
  res <- confidence_interval(X,c(200,300))
  expect_equal(res$change,c(200, 300))
  expect_equal(round(res$lower,4),c(37.5143, 237.1052))
  expect_equal(round(res$upper,4),c(382.2948, 370.5625))

  expect_error(confidence_interval(X,c()))
  expect_error(confidence_interval(X,c(0)))
})
