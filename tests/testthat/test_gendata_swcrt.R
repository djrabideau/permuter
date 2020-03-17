context("Testing gendata_swcrt")

test_that("gendata_swcrt has correct number of rows", {
  ds <- gendata_swcrt(nclus = 6, size = 10, nstep = 2, theta = 0,
                      sigma = 0.5, mu = 0, beta = 0:2, nu = 0.2, sd = 1)
  expect_equal(nrow(ds), 180)
})

test_that("gendata_swcrt errors", {
  expect_error( # nclus not a multiple of nstep
    gendata_swcrt(nclus = 6, size = 10, nstep = 4, theta = 0,
                      sigma = 0.5, mu = 0, beta = 0:2, nu = 0.2, sd = 1)
  )
  expect_error( # beta wrong length
    gendata_swcrt(nclus = 6, size = 10, nstep = 2, theta = 0,
                      sigma = 0.5, mu = 0, beta = 0:3, nu = 0.2, sd = 1)
  )
})
