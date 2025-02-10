test_that("Symmetric Confirmation", {
  acv <- autocovariance(electricity,0)
  aco <- autocorrelation(electricity,0)

  expect_equal(acv, t(acv))
  expect_equal(aco, t(aco))
  expect_equal(round(diag(aco),7),
               c(0.8847832, 0.8831488, 0.8893171, 0.8978068, 0.8821414,
                 0.8605953, 0.9524064, 1.0406137, 1.1376170, 1.1286609,
                 1.0840598, 1.0939983, 1.0963860, 1.1010148, 1.0733718,
                 1.0808826, 1.1329568, 1.1108067, 1.0396160, 1.0130921,
                 0.9859969, 0.8569262, 0.8123149, 0.8077560))
})

test_that("Some autocorrelation on electricity", {
  aco <- autocorrelation(electricity,0)
  aco1 <- autocorrelation(electricity,1)

  expect_equal(round(diag(aco),7),
               c(0.8847832, 0.8831488, 0.8893171, 0.8978068, 0.8821414,
                 0.8605953, 0.9524064, 1.0406137, 1.1376170, 1.1286609,
                 1.0840598, 1.0939983, 1.0963860, 1.1010148, 1.0733718,
                 1.0808826, 1.1329568, 1.1108067, 1.0396160, 1.0130921,
                 0.9859969, 0.8569262, 0.8123149, 0.8077560))
  expect_equal(round(diag(aco1),7),
               c(0.6283851, 0.6354570, 0.6416834, 0.6605441, 0.6517430,
                 0.6430812, 0.7043002, 0.6972579, 0.7445781, 0.6944436,
                 0.7423818, 0.7947165, 0.7950447, 0.8375397, 0.8159054,
                 0.8189858, 0.8597472, 0.8248445, 0.7173009, 0.5845091,
                 0.5641858, 0.5248945, 0.6026876, 0.6147604))
})

test_that("Lengths and Sizes", {
  acv <- autocovariance(electricity,0)
  acv1 <- autocovariance(electricity,0:1)
  aco <- autocorrelation(electricity,0)
  aco1 <- autocorrelation(electricity,0:1)

  expect_equal(dim(acv),c(nrow(electricity),nrow(electricity)))
  expect_equal(length(acv),nrow(electricity)^2)

  expect_equal(dim(acv1),NULL)
  expect_equal(length(acv1),2)


  expect_equal(dim(aco),c(nrow(electricity),nrow(electricity)))
  expect_equal(length(aco),nrow(electricity)^2)

  expect_equal(dim(aco1),NULL)
  expect_equal(length(aco1),2)
})
