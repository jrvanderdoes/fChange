## SINGLE CHANGES
test_that("Characteristic change", {
  if(Sys.info()['sysname'] =='Linux') {
    expect_false(FALSE)
    return()
  }
  ## No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- fchange(X = dat, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.835)

  ## Change
  set.seed(123)
  tmp <- fchange(electricity, method='characteristic', h=0,
                 statistic='Mn', type='Single',
                 critical = 'simulation',
                 M = 200)
  expect_equal(tmp$pvalue, 0)


  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################
  ## No Change
  set.seed(123)
  tmp <- fchange(dat, method='characteristic', h=0,
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.755)

  set.seed(123)
  tmp <- fchange(X = dat, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'welch',
                M = 200)
  expect_equal(tmp$pvalue, 0.810)
  expect_error(change(dat, method='characteristic', h=0,
                      statistic='Mn', type='Single',
                      critical = 'welch',
                      M = 200))

  set.seed(123)
  tmp <- fchange(X = dat, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'resample',
                M = 200)
  expect_equal(tmp$pvalue, 0.880)
  tmp <- fchange(dat, method='characteristic', h=0,
                statistic='Mn', type='Single',
                critical = 'resample',
                M = 200)
  expect_equal(tmp$pvalue, 0.700)


  ## Change
  set.seed(123)
  tmp <- fchange(X = electricity, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0)
  expect_equal(tmp$pvalue, 0)

  set.seed(123)
  tmp <- fchange(X = electricity, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'welch',
                M = 200)
  # expect_equal(tmp$pvalue, 0)
  expect_error(change(electricity, method='characteristic', h=0,
                      statistic='Mn', type='Single',
                      critical = 'welch',
                      M = 200))

  set.seed(123)
  tmp <- fchange(X = electricity, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'resample',
                M = 200)
  expect_equal(tmp$pvalue, 0)
  tmp <- fchange(electricity, method='characteristic', h=0,
                statistic='Mn', type='Single',
                critical = 'resample',
                M = 200)
  expect_equal(tmp$pvalue, 0)
})


test_that("Mean change", {
  # Simulation - No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- fchange(dat, method='mean',
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0.511)

  tmp <- fchange(electricity, method='mean',
                 statistic='Mn', type='Single',
                 critical = 'resample',
                 M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0)

  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################

  # Simulation - No Change
  set.seed(123)
  tmp <- fchange(dat, method='mean',
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0.579)

  # resample - Change
  set.seed(1234)
  tmp <- fchange(electricity, method='mean',
                statistic='Tn', type='Single',
                critical = 'resample',
                M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0)
})


test_that("Robust change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- fchange(dat, method='robustmean',
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.493)

  # Change
  set.seed(134)
  dat_change <- dat
  dat_change$data[,50:100] <- dat_change$data[,50:100]+5

  tmp <- fchange(dat_change, method='robustmean',
                 statistic='Mn', type='Single',
                 critical = 'resample',
                 M = 500)
  expect_equal(tmp$pvalue, 0.014)

  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################
  # No Change
  set.seed(123)
  tmp <- fchange(dat, method='robustmean',
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.778)

  # Change
  set.seed(213)
  dat_change <- dat
  dat_change$data[,50:100] <- dat_change$data[,50:100]+5

  tmp <- fchange(dat_change, method='robustmean',
                statistic='Tn', type='Single',
                critical = 'resample',
                M = 500)
  expect_equal(tmp$pvalue, 0.008)
})


test_that("Eigen change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- fchange(dat, method='eigenjoint', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.240)

  tmp <- fchange(electricity, method='eigensingle', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)

  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################
  ## No Change
  set.seed(123)
  tmp <- fchange(dat, method='eigensingle', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.210)

  set.seed(123)
  tmp <- fchange(dat, method='eigenjoint', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'resample',
                M = 5)
  expect_equal(tmp$pvalue, 0.4)
  tmp <- fchange(dat, method='eigensingle', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'resample',
                M = 5)
  expect_equal(tmp$pvalue, 0)

  ## Change
  set.seed(123)
  tmp <- fchange(electricity, method='eigenjoint', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)
  tmp <- fchange(electricity, method='eigensingle', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)

  tmp <- fchange(electricity, method='eigenjoint', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'resample',
                M = 2)
  expect_equal(tmp$pvalue, 0)
  tmp <- fchange(electricity, method='eigensingle', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'resample',
                M = 10)
  expect_equal(tmp$pvalue,0)
})


test_that("Trace change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- fchange(dat, method='trace',
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.292)


  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################
  # No change
  tmp <- fchange(dat, method='trace',
                statistic='Mn', type='Single',
                critical = 'resample',
                M = 1000)
  expect_equal(tmp$pvalue, 0.150)

  # Change
  tmp <- fchange(electricity, method='trace',
                statistic='Tn', type='Single',
                critical = 'resample',
                M = 1000)
  expect_equal(tmp$pvalue, 0.016)

  tmp <- fchange(electricity, method='trace',
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.011)
})


test_that("Covariance change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- fchange(dat, method='covariance',
                statistic='Tn',
                critical = 'simulation',
                M = 1000, cov.res = 20)
  expect_equal(tmp$pvalue, 0.051)


  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################

  set.seed(123)
  tmp <- fchange(dat, method='covariance',
                statistic='Tn',
                critical = 'resample',
                M = 10, cov.res = 5)
  expect_equal(tmp$pvalue, 0.2)

  expect_error(change(dat, method='covariance',
                statistic='Mn',
                critical = 'simulation',
                M = 1000))

  # Change
  set.seed(123)
  tmp <- fchange(electricity, method='covariance',
                statistic='Tn',
                critical = 'simulation',
                M = 1000, cov.res = 15)
  expect_equal(tmp$pvalue, 0)
})


test_that("PCA-Mean change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- fchange(dat, method='projmean',
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.210)


  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################

  # No Change
  tmp <- fchange(dat, method='projmean',
                statistic='Tn',
                critical = 'resample',
                M = 1000)
  expect_equal(tmp$pvalue, 0.067)

  tmp <- fchange(dat, method='projmean',
                statistic='Mn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.244)
  tmp <- fchange(dat, method='projmean',
                statistic='Mn',
                critical = 'resample',
                M = 1000)
  expect_equal(tmp$pvalue, 0.135)

  # Change
  set.seed(123)
  tmp <- fchange(electricity, method='projmean',
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)
  tmp <- fchange(electricity, method='projmean',
                statistic='Mn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)
})


test_that("PCA-Distribution change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(50)
  tmp <- fchange(dat, method='projdistribution',
                statistic='Tn',
                critical = 'resample',
                M = 100)
  expect_equal(tmp$pvalue, 0.530)
  tmp <- fchange(dat, method='projdistribution',
                 statistic='Mn',
                 critical = 'resample',
                 M = 100)
  expect_equal(tmp$pvalue, 0.53)


  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################

  # No Change
  set.seed(123)
  tmp <- fchange(dat, method='projdistribution',
                statistic='Mn',
                critical = 'resample',
                M = 100)
  expect_equal(tmp$pvalue, 0.610)

  # Change
  dat$data[,25:50] <- dat$data[,25:50]+5
  tmp <- fchange(dat, method='projdistribution',
                statistic='Tn',
                critical = 'resample',
                M = 100)
  expect_equal(tmp$pvalue, 0)


  # Error
  expect_error(change(dat, method='projdistribution',
                      statistic='Mn',
                      critical = 'simulation',
                      M = 100))
})


## ELBOW PLOTS
test_that("Check Elbow Plots", {
  ## No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  res <- fchange(X = dat,
                method='characteristic',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)
  expect_equal(class(res$suggestion$plot)[1], 'gg')
  expect_equal(class(res$plots$variability)[1], 'gg')
  expect_equal(class(res$plots$explained)[1], 'gg')
  expect_equal(class(res$plots$improvement)[1], 'gg')
  expect_equal(round(sum(res$information[,2]),4), 35.6545)

  ## Change
  set.seed(123)
  dat_change <- generate_brownian_bridge(100)
  dat_change$data[,51:100] <- dat_change$data[,51:100] + 5
  res <- fchange(X = dat_change,
                 method='characteristic',
                 type='Elbow',
                 max_changes=min(ncol(dat),20),
                 eigen_number=2, h=1,
                 W = space_measuring_functions(X = dat, M = 20, space='BM'),
                 K = bartlett_kernel,
                 weighting = 1/4, TVE=0.95,
                 trim_function = function(X) { 0 },
                 errors='L2', recommendation_change_points = 2,
                 recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)
  expect_equal(class(res$suggestion$plot)[1], 'gg')
  expect_equal(class(res$plots$variability)[1], 'gg')
  expect_equal(class(res$plots$explained)[1], 'gg')
  expect_equal(class(res$plots$improvement)[1], 'gg')
  expect_equal(round(sum(res$information[,2]),4),  224.252)

  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################
  ## No Changes

  set.seed(123)
  res <- fchange(X = dat,
                method='mean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat,
                method='robustmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat,
                method='eigenjoint',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat,
                method='eigensingle',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat,
                method='trace',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat,
                method='covariance',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat,
                method='projmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat,
                method='projdistribution',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)



  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################

  ## Changes
  set.seed(123)
  res <- fchange(X = dat_change,
                method='mean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  res <- fchange(X = dat_change,
                method='robustmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  res <- fchange(X = dat_change,
                method='eigenjoint',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,c(52, 49, 50))

  res <- fchange(X = dat_change,
                method='eigensingle',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,c(52, 49, 50))

  res <- fchange(X = dat_change,
                method='trace',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- fchange(X = dat_change,
                method='covariance',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,c(39,50))

  res <- fchange(X = dat_change,
                method='projmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  res <- fchange(X = dat_change,
                method='projdistribution',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  ## Errors
  set.seed(1234)
  dat <- generate_brownian_bridge(100)
  res <- fchange(X = dat,
                method='characteristic',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='trace', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  dat$data[,51:100] <- dat$data[,51:100] +5
  res <- fchange(X = dat,
                method='characteristic',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_functions(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='trace', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)
})


## Binary Segmentation
test_that("Binary Segmentation", {
  set.seed(12345)
  X <- generate_brownian_bridge(150)
  X[,51:150] <- X[,51:150] + 3
  X[,101:150] <- X[,101:150] * 2

  res <- fchange(X,
                method='characteristic',
                statistic='Tn',
                critical='welch',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,c(50,99))
  expect_equal(res$pvalue,c(0.000,0.000))

  #################
  ## Only run below locally
  testthat::skip_on_cran()
  #################

  res <- fchange(X,
                method='mean',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,c(50,66,100))
  expect_equal(res$pvalue,c(0.001,0,0))


  res <- fchange(X,
                method='robustmean',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,c(50,66,100))
  expect_equal(res$pvalue,c(0,0,0))

  res <- fchange(X,
                method='eigenjoint',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,c(102))
  expect_equal(res$pvalue,c(0))

  res <- fchange(X,
                method='eigensingle',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,c(102))
  expect_equal(res$pvalue,c(0.003))

  res <- fchange(X,
                method='trace',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,NULL)
  expect_equal(res$pvalue,NULL)

  res <- fchange(X,
                method='projmean',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,c(50,66,100))
  expect_equal(res$pvalue,c(0,0,0))


  set.seed(123)
  X <- generate_brownian_bridge(50)
  X$data[,26:50] <- X$data[,26:50] +5

  res <- fchange(X,
                method='projdistribution',
                statistic='Tn',
                critical='resample',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,c(25))
  expect_equal(res$pvalue,c(0))

  res <- fchange(X,
                method='covariance',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                resample_blocks = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_functions(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                silent.binary = TRUE)
  expect_equal(res$location,NULL)
  expect_equal(res$pvalue,NULL)
})
