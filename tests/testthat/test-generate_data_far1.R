test_that("FAR(1) Test", {
  set.seed(1234)
  dat <- generate_far1(N=10, resolution=4, sd=1, dependence=0, drop_first=TRUE)
  expect_equal(dim(dat$data), c(3,10))
  expect_equal(round(dat$data[1,2:4],4), c(0.1602, 0.6261, -1.3543))


  dat <- generate_far1(N=10, resolution=c(0,0.2,0.5,1),
                       sd=1, dependence=0, drop_first=FALSE)
  expect_equal(dim(dat$data), c(4,10))
  expect_equal(sum(dat$data[1,]), 0)
  expect_equal(round(dat$data[2,4],4), -0.2242)
})
