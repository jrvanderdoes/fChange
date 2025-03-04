test_that("Brownian Motion", {
  set.seed(234)
  tmp <- generate_brownian_motion(10,v=20,sd=1)
  expect_equal(dim(tmp), c(20,10))
  expect_equal(sum(tmp$data[1,]), 0)

  set.seed(234)
  tmp <- generate_brownian_motion(10,v=c(0,0.25,0.75,1),sd=1)
  expect_equal(round(tmp$data[2,4],4),0.7356)
})

test_that("Brownian Bridge",{
  set.seed(234)
  tmp <- generate_brownian_bridge(10,v=20,sd=2)
  expect_equal(dim(tmp), c(20,10))
  expect_equal(sum(tmp$data[1,]), 0)

  set.seed(234)
  tmp <- generate_brownian_bridge(10,v=c(0,0.25,0.75,1),sd=1)
  expect_equal(round(tmp$data[2,4],4),0.5489)
})
