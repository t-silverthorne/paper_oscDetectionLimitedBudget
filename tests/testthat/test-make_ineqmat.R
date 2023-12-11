test_that("warning generated", {
  # should warn that ineqmat is not to be used for convex programming
  opts = make_default_opts(prob_size='small')
  expect_warning(make_ineqmat(opts),
                 'The output of this function is only useful for enforcing lattice constraints.*')
})
