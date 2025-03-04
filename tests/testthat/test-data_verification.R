test_that('Verify Inputs', {
  expect_equal(.verify_input('A',c('A','B','C')), 'A')
  expect_equal(.verify_input('a',c('A','B','C')), 'A')
  expect_equal(.verify_input('A',c('a','b','c')), 'a')

  expect_error(.verify_input('A',c('B','C')))

  expect_equal(.verify_input(c('B','A'),c('A','B','C')), 'B')
  expect_equal(.verify_input(c('D','A'),c('A','B','C')), 'A')
})
