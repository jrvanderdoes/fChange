test_that("CIDR and OCIDR", {
  tmp <- dfts(SPYUS500$data[,1:100], name='SP500 100 Days',
   labels = SPYUS500$labels[1:100], fparam=SPYUS500$fparam)
  c_vals <- cidr(tmp)
  o_vals <- ocidr(tmp)
  expect_equal(round(sum(c_vals$data),3), 2785.417)
  expect_equal(sum(c_vals$data[1,]^2), 0)


  expect_equal(round(sum(o_vals$data),3), 3005.710)
  expect_equal(round(sum(o_vals$data[1,]^2),3), 33.369)
})
