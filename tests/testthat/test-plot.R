test_that("Plot Options", {
  tmp <- plot.dfts(electricity, changes=c(50,175),type='spaghetti')
  expect_equal(class(tmp)[1],"plotly")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='highdim')
  expect_equal(class(tmp)[1],"trellis")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='rainbow')
  expect_equal(class(tmp)[1],"gg")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='banded')
  expect_equal(class(tmp)[1],"gg")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='acf')
  expect_equal(class(tmp)[1],"list")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='pacf')
  expect_equal(class(tmp)[1],"list")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='summary')
  expect_equal(class(tmp)[1],"list")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='qq')
  expect_equal(class(tmp)[1],"gg")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='distribution')
  expect_equal(class(tmp)[1],"gg")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='change')
  expect_equal(class(tmp)[1],"plotly")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='interval',int.gradual = F)
  expect_equal(class(tmp)[1],"plotly")
  tmp <- plot.dfts(electricity, changes=c(50,175),type='interval',int.gradual = T)
  expect_equal(class(tmp)[1],"plotly")

  tmp <- plot.dfts(electricity, changes=c(50,175),type='surface')
  expect_equal(class(tmp)[1],"matrix")
})
