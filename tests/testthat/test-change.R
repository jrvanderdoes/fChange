## SINGLE CHANGES
test_that("Characteristic change", {
  ## No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- change(X = dat, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.960)
  tmp <- change(dat, method='characteristic', h=0,
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.995)

  set.seed(123)
  tmp <- change(X = dat, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'welch',
                M = 200)
  expect_equal(tmp$pvalue, 0.870)
  expect_error(change(dat, method='characteristic', h=0,
                      statistic='Mn', type='Single',
                      critical = 'welch',
                      M = 200))

  set.seed(123)
  tmp <- change(X = dat, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'permutation',
                M = 200)
  expect_equal(tmp$pvalue, 0.305)
  tmp <- change(dat, method='characteristic', h=0,
                statistic='Mn', type='Single',
                critical = 'permutation',
                M = 200)
  expect_equal(tmp$pvalue, 0.520)


  ## Change
  set.seed(123)
  tmp <- change(X = electricity, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0)
  tmp <- change(electricity, method='characteristic', h=0,
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0)

  set.seed(123)
  tmp <- change(X = electricity, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'welch',
                M = 200)
  expect_equal(tmp$pvalue, 0)
  expect_error(change(electricity, method='characteristic', h=0,
                      statistic='Mn', type='Single',
                      critical = 'welch',
                      M = 200))

  set.seed(123)
  tmp <- change(X = electricity, method='characteristic', h=0,
                statistic='Tn', type='Single',
                critical = 'permutation',
                M = 200)
  expect_equal(tmp$pvalue, 0)
  tmp <- change(electricity, method='characteristic', h=0,
                statistic='Mn', type='Single',
                critical = 'permutation',
                M = 200)
  expect_equal(tmp$pvalue, 0)
})


test_that("Mean change", {
  # Simulation - No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- change(dat, method='mean',
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0.369)

  tmp <- change(dat, method='mean',
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0.571)

  # Permutation - Change
  tmp <- change(electricity, method='mean',
                statistic='Tn', type='Single',
                critical = 'permutation',
                M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0)

  tmp <- change(electricity, method='mean',
                statistic='Mn', type='Single',
                critical = 'permutation',
                M = 1000, K = bartlett_kernel)
  expect_equal(tmp$pvalue, 0)
})


test_that("Robust change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- change(dat, method='robustmean',
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.263)

  tmp <- change(dat, method='robustmean',
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.418)

  # Change
  dat$data[,50:100] <- dat$data[,50:100]+5
  tmp <- change(dat, method='robustmean',
                statistic='Tn', type='Single',
                critical = 'permutation',
                M = 500)
  expect_equal(tmp$pvalue, 0.012)

  tmp <- change(dat, method='robustmean',
                statistic='Mn', type='Single',
                critical = 'permutation',
                M = 500)
  expect_equal(tmp$pvalue, 0.012)
})


test_that("Eigen change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- change(dat, method='eigenjoint', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.615)
  tmp <- change(dat, method='eigensingle', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 200)
  expect_equal(tmp$pvalue, 0.63)

  set.seed(123)
  # dat <- generate_brownian_bridge(30)
  tmp <- change(dat, method='eigenjoint', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'permutation',
                M = 5)
  expect_equal(tmp$pvalue, 0.6)
  tmp <- change(dat, method='eigensingle', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'permutation',
                M = 5)
  expect_equal(tmp$pvalue, 0.4)

  ## Change
  set.seed(123)
  tmp <- change(electricity, method='eigenjoint', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)
  tmp <- change(electricity, method='eigensingle', h=0, eigen_number=1,
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)

  tmp <- change(electricity, method='eigenjoint', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'permutation',
                M = 2)
  expect_equal(tmp$pvalue, 0)
  tmp <- change(electricity, method='eigensingle', h=0, eigen_number=1,
                statistic='Mn',
                critical = 'permutation',
                M = 10)
  expect_equal(tmp$pvalue,0)
})


test_that("Trace change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- change(dat, method='trace',
                statistic='Tn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.487)

  tmp <- change(dat, method='trace',
                statistic='Mn', type='Single',
                critical = 'permutation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.401)

  # Change
  tmp <- change(electricity, method='trace',
                statistic='Tn', type='Single',
                critical = 'permutation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.01)

  tmp <- change(electricity, method='trace',
                statistic='Mn', type='Single',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.009)
})


test_that("Covariance change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- change(dat, method='covariance',
                statistic='Tn',
                critical = 'simulation',
                M = 1000, cov.res = 20)
  expect_equal(tmp$pvalue, 0.651)
  set.seed(123)
  tmp <- change(dat, method='covariance',
                statistic='Tn',
                critical = 'permutation',
                M = 10, cov.res = 5)
  expect_equal(tmp$pvalue, 0.6)

  expect_error(change(dat, method='covariance',
                statistic='Mn',
                critical = 'simulation',
                M = 1000))

  # Change
  set.seed(123)
  tmp <- change(electricity, method='covariance',
                statistic='Tn',
                critical = 'simulation',
                M = 1000, cov.res = 15)
  expect_equal(tmp$pvalue, 0)
})


test_that("PCA-Mean change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(100)
  tmp <- change(dat, method='projmean',
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.327)
  tmp <- change(dat, method='projmean',
                statistic='Tn',
                critical = 'permutation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.182)

  tmp <- change(dat, method='projmean',
                statistic='Mn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.338)
  tmp <- change(dat, method='projmean',
                statistic='Mn',
                critical = 'permutation',
                M = 1000)
  expect_equal(tmp$pvalue, 0.208)

  # Change
  set.seed(123)
  tmp <- change(electricity, method='projmean',
                statistic='Tn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)
  tmp <- change(electricity, method='projmean',
                statistic='Mn',
                critical = 'simulation',
                M = 1000)
  expect_equal(tmp$pvalue, 0)
})


test_that("PCA-Distribution change", {
  # No Change
  set.seed(123)
  dat <- generate_brownian_bridge(50)
  tmp <- change(dat, method='projdistribution',
                statistic='Tn',
                critical = 'permutation',
                M = 100)
  expect_equal(tmp$pvalue, 0.630)
  tmp <- change(dat, method='projdistribution',
                statistic='Mn',
                critical = 'permutation',
                M = 100)
  expect_equal(tmp$pvalue, 0.870)

  # Change
  dat$data[,25:50] <- dat$data[,25:50]+5
  tmp <- change(dat, method='projdistribution',
                statistic='Tn',
                critical = 'permutation',
                M = 100)
  expect_equal(tmp$pvalue, 0)
  tmp <- change(dat, method='projdistribution',
                statistic='Mn',
                critical = 'permutation',
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
  res <- change(X = dat,
                method='characteristic',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
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
  expect_equal(round(sum(res$information[,2]),4), 178.2479)

  res <- change(X = dat,
                method='mean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='robustmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='eigenjoint',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='eigensingle',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='trace',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='covariance',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='projmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='projdistribution',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)


  ## Change
  set.seed(123)
  dat$data[,51:100] <- dat$data[,51:100] + 5
  res <- change(X = dat,
                method='characteristic',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
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
  expect_equal(round(sum(res$information[,2]),4),  820.6515)

  res <- change(X = dat,
                method='mean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  res <- change(X = dat,
                method='robustmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  res <- change(X = dat,
                method='eigenjoint',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,c(9,74,34,55,50))

  res <- change(X = dat,
                method='eigensingle',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,c(9,74,34,68,59,54,49,52,50))

  res <- change(X = dat,
                method='trace',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  res <- change(X = dat,
                method='covariance',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,c(92,49,50))

  res <- change(X = dat,
                method='projmean',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  res <- change(X = dat,
                method='projdistribution',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='L2', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)

  ## Errors
  set.seed(1234)
  dat <- generate_brownian_bridge(100)
  res <- change(X = dat,
                method='characteristic',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='trace', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,NA)

  dat$data[,51:100] <- dat$data[,51:100] +5
  res <- change(X = dat,
                method='characteristic',
                type='Elbow',
                max_changes=min(ncol(dat),20),
                eigen_number=2, h=1,
                W = space_measuring_vectors(X = dat, M = 20, space='BM'),
                K = bartlett_kernel,
                weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 },
                errors='trace', recommendation_change_points = 2,
                recommendation_improvement = 0.15)
  expect_equal(res$suggestion$changes,50)
})


## Binary Segmentation
test_that("Binary Segmentation", {
  set.seed(123)
  X <- generate_brownian_bridge(150)
  X$data[,51:150] <- X$data[,51:150] +5
  X$data[,101:150] <- X$data[,101:150]*2


  res <- change(X,
                method='characteristic',
                statistic='Tn',
                critical='welch',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,c(49,99))
  expect_equal(res$pvalue,c(0.034,0.027))

  res <- change(X,
                method='mean',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,c(50,100))
  expect_equal(res$pvalue,c(0.006,0))

  res <- change(X,
                method='robustmean',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,c(50,100))
  expect_equal(res$pvalue,c(0,0.027))

  res <- change(X,
                method='eigenjoint',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,c(101))
  expect_equal(res$pvalue,c(0))

  res <- change(X,
                method='eigensingle',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,c(99))
  expect_equal(res$pvalue,c(0))

  res <- change(X,
                method='trace',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,NULL)
  expect_equal(res$pvalue,NULL)

  res <- change(X,
                method='projmean',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,c(50,100))
  expect_equal(res$pvalue,c(0,0))


  set.seed(123)
  X <- generate_brownian_bridge(50)
  X$data[,26:50] <- X$data[,26:50] +5

  res <- change(X,
                method='projdistribution',
                statistic='Tn',
                critical='permutation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,c(25))
  expect_equal(res$pvalue,c(0))

  res <- change(X,
                method='covariance',
                statistic='Tn',
                critical='simulation',
                type='segmentation',
                perm_type = 'separate', replace=TRUE,
                blocksize = 1,
                eigen_number=3, h=3,
                M = 1000, J=50,
                W = space_measuring_vectors(X = X, M = 20, space='BM'),
                K = bartlett_kernel,
                alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                trim_function = function(X) { 0 })
  expect_equal(res$location,NULL)
  expect_equal(res$pvalue,NULL)
})
