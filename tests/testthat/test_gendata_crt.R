context("Testing gendata_crt")

test_that("gendata_crt has correct dimensions", {
  ds <- gendata_crt(nclus = c(10, 10), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  expect_equal(dim(ds), c(200, 5))
})

test_that("gendata_crt errors", {
  expect_error( # nclus not vector of 2
    gendata_crt(nclus = 10, size = c(10, 10), sigma = 1, mu = 0, sd = 1)
  )
  expect_error( # size not a vector of 2
    gendata_crt(nclus = c(10, 10), size = 10, sigma = 1, mu = 0, sd = 1)
  )
  expect_error(
    gendata_crt(family = Gamma, nclus = c(2, 2), size = c(2, 2), sigma = 1,
                mu = 0)
  )
})
