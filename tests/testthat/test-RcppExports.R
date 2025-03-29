test_that("Integrate Verify", {
  expect_equal(dot_integrate(1:10), 5.5)
  expect_equal(
    dot_integrate(1:10, c(0,0.1,seq(0.5,1,length.out=8))), 4.4)
})

test_that("Integrate Columns Verify", {
  dat <- matrix(c(1:10,c(rep(5,5),rep(10,5))),ncol = 2,nrow = 10,byrow = FALSE)

  expect_equal( dot_integrate_col(dat), c(5.5,7.5))
  expect_equal(
    round(dot_integrate_col(dat, c(0,0.1,seq(0.5,1,length.out=8))),7),
    c(4.400000, 6.6071429))
})

test_that("Integrate Agreement Verify", {
  dat <- matrix(c(1:10,c(rep(2,5),rep(8,5))),ncol = 2,nrow = 10,byrow = FALSE)
  ts <- c(seq(0,0.5, length.out=7),seq(0.7,1,length.out=3))

  expect_equal( dot_integrate_col(dat),
                c(dot_integrate(dat[,1]),dot_integrate(dat[,2])))
  expect_equal( dot_integrate_col(dat,ts),
                c(dot_integrate(dat[,1],ts),dot_integrate(dat[,2],ts)))
})
