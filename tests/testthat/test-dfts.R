test_that("Verify Input Formats", {
  dfts(dfts(electricity))
  if(requireNamespace("fda", quietly = TRUE)){
    tmp <- dfts(fda::Data2fd(1:24,electricity$data))
    expect_equal(tmp$data[2,7],electricity$data[2,7])
  }
  if(requireNamespace("fda.usc", quietly = TRUE)){
    tmp <- dfts(fda.usc::fdata(t(electricity$data)))
    expect_equal(tmp$data[2,7],electricity$data[2,7])
  }
  if(requireNamespace("rainbow", quietly = TRUE)){
    tmp <- dfts(suppressWarnings(rainbow::fts(x=1:nrow(electricity$data),y=electricity$data)))
    expect_equal(tmp$data[2,7],electricity$data[2,7])
    tmp <- dfts(rainbow::fds(1:nrow(electricity$data),electricity$data))
    expect_equal(tmp$data[2,7],electricity$data[2,7])
  }
  if(requireNamespace('funData', quietly = TRUE) ){
    tmp <- dfts( funData::funData(1:24, t(electricity$data)) )
    expect_equal(tmp$data[2,7],electricity$data[2,7])
  }

  tmp <- dfts(electricity)
  expect_equal(tmp$data[2,7],electricity$data[2,7])
})

test_that("Check inheritance of dfts", {
  expect_false(is.dfts(electricity$data))
  expect_true(is.dfts(electricity))
})

test_that("Check math of dfts", {
  elec <- dfts(electricity)

  expect_equal((elec+elec)$data[8,7],electricity$data[8,7]*2)
  expect_equal(sum((elec-elec)$data),0)
  expect_equal(round(sqrt(elec)$data[8,9],5),6.64605)

  expect_equal(min(elec)$data[1:3], c(0.2,0.0,0.0) )
  expect_equal(max(elec)$data[1:3], c(67.2,64.6,62.68) )

  expect_equal(mean(elec)[15], mean(electricity$data[15,]) )
  expect_equal(median(elec)[8], median(electricity$data[8,]) )
  expect_equal(quantile(elec,probs = 0.95)[1,],
               quantile(electricity$data[1,],prob=0.95)[[1]])
})

test_that("Check lag of dfts", {
  elec <- dfts(electricity)

  expect_equal(dim(lag(elec,1)),c(24,364))
  expect_equal(lag(elec,lag=1,difference=2)$data[22,42],18.69)
})
