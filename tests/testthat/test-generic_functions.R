test_that("Center works", {
  expect_equal(center(1:3), c(-1,0,1))
  expect_equal(matrix(center(data.frame(rep(1,3),rep(2,3),rep(3,3))),3,3),
               matrix(rep(-1:1,time2=3),3,3,byrow = TRUE))
  expect_equal(matrix(center(dfts(data.frame(rep(1,3),rep(2,3),rep(3,3))))$data,3,3),
               matrix(rep(-1:1,time2=3),3,3,byrow = TRUE))
})

test_that("PCA works", {
  expect_equal(round(pca(1:10)$sdev,5), 3.02765)
  expect_equal(round(pca(electricity)$sdev[1:3],6),
               c(15.764828, 5.456613, 3.313120))
})

test_that("SD works", {
  expect_equal(round(sd(1:10),5), 3.02765)
  expect_equal(round(sd(electricity)[1:3],5),
               c(16.35713, 16.34201, 16.39898))
})

test_that("Var works", {
  expect_equal(round(var(1:10),5), 9.16667)
  expect_equal(round(var(electricity,type='pointwise')[1:3],5),
               c(267.55565, 267.06142, 268.92667))
  expect_equal(round(diag(var(electricity,type='operator'))[1:3],5),
               c(266.82262, 266.32974, 268.18988))
})

test_that("CUSUM works", {
  expect_equal(round(cumsum(electricity)$data[c(1,5,10),365],5),
               c(14562.59, 10548.78, 16999.70))
})










