context("Testing permtest")

test_that("permtest seed parameter works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)

  # not parallel
  tmp1 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = 1, seed = 123, quietly = T)
  tmp2 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = 1, seed = 123, quietly = T)
  tmp3 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = 1, quietly = T)
  tmp4 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = 1, seed = 124, quietly = T)
  expect_true(identical(tmp1, tmp2))
  expect_false(identical(tmp2, tmp3))
  expect_false(identical(tmp2, tmp4))

  # parallel
  ncores <- parallel::detectCores() - 1
  tmp1 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = ncores, seed = 123, quietly = T)
  tmp2 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = ncores, seed = 123, quietly = T)
  tmp3 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = ncores, quietly = T)
  tmp4 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = ncores, seed = 124, quietly = T)
  expect_true(identical(tmp1, tmp2))
  expect_false(identical(tmp2, tmp3))
  expect_false(identical(tmp2, tmp4))
})

test_that("permtest with user defined function works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  # permtest_glm
  tmp5 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 99, ncores = 1, seed = 444, quietly = T)
  # permtest
  f_tmp <- function(d) {
    fit_tmp <- glm(y ~ trt, data = d)
    coef_tmp <- coef(fit_tmp)['trt']
    return(as.numeric(coef_tmp))
  }
  tmp6 <- permtest(f_tmp, trtname = 'trt', runit = 'clusid', data = ds,
                   nperm = 99, ncores = 1, seed = 444, quietly = T)
  expect_true(identical(tmp5, tmp6))
})

# should add other permtest_xxx tests, but first need dummy data
