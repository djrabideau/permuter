context("Testing permtest")

test_that("permtest parallel works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  tmp1 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, quietly = T)
  expect_true(getDoParWorkers() == 1)
  tmp2 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 2, quietly = T)
  expect_true(getDoParWorkers() == 2)
  tmp3 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, quietly = T)
  expect_true(getDoParWorkers() == 1)
})

test_that("permtest seed parameter works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)

  # not parallel
  tmp1 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, seed = 123, quietly = T)
  tmp2 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, seed = 123, quietly = T)
  tmp3 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, quietly = T)
  tmp4 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, seed = 124, quietly = T)
  expect_true(identical(tmp1$permCoefs, tmp2$permCoefs))
  expect_false(identical(tmp2$permCoefs, tmp3$permCoefs))
  expect_false(identical(tmp2$permCoefs, tmp4$permCoefs))

  # parallel
  ncores <- parallel::detectCores() - 1
  tmp1 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = ncores, seed = 123, quietly = T)
  tmp2 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = ncores, seed = 123, quietly = T)
  tmp3 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = ncores, quietly = T)
  tmp4 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = ncores, seed = 124, quietly = T)
  expect_true(identical(tmp1$permCoefs, tmp2$permCoefs))
  expect_false(identical(tmp2$permCoefs, tmp3$permCoefs))
  expect_false(identical(tmp2$permCoefs, tmp4$permCoefs))

  # not parallel vs parallel
  tmp1 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, seed = 123, quietly = T)
  tmp2 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 2, seed = 123, quietly = T)
  expect_true(identical(tmp1$permCoefs, tmp2$permCoefs))
})

test_that("permtest with fitted model object works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  # permtest_glm
  tmp5 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, seed = 444, quietly = T)
  # permtest
  m1 <- glm(y ~ trt, data = ds)
  tmp6 <- permtest(m1, trtname = 'trt', runit = 'clusid', data = ds,
                   nperm = 10, ncores = 1, seed = 444, quietly = T)
  expect_true(identical(tmp5$permCoefs, tmp6$permCoefs))
})

test_that("permtest with user defined function works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  # permtest_glm
  tmp7 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 10, ncores = 1, seed = 444, quietly = T)
  # permtest
  f_tmp <- function(d) {
    fit_tmp <- glm(y ~ trt, data = d)
    coef_tmp <- coef(fit_tmp)['trt']
    return(as.numeric(coef_tmp))
  }
  tmp8 <- permtest(f_tmp, trtname = 'trt', runit = 'clusid', data = ds,
                   nperm = 10, ncores = 1, seed = 444, quietly = T)
  expect_true(identical(tmp7$permCoefs, tmp8$permCoefs))
})

test_that("permtest with non-zero null works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  # permtest
  m1 <- glm(y ~ trt, data = ds)
  expect_no_condition(  tmp9 <- permtest(m1, trtname = 'trt', runit = 'clusid', data = ds,
                                         alternative = 'less', theta = 1,
                                         nperm = 10, ncores = 1, seed = 444, quietly = T))
  expect_error(  tmp10 <- permtest({function(x) x}, trtname = 'trt', runit = 'clusid', data = ds,
                                         alternative = 'less', theta = 1,
                                         nperm = 10, ncores = 1, seed = 444, quietly = T))
})

test_that("permtest_glm with non-zero null works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  # permtest
  expect_no_condition(  tmp11 <- permtest_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                                         alternative = 'less', theta = 1,
                                         nperm = 10, ncores = 1, seed = 444, quietly = T))
})

# should add other permtest_xxx tests, but first need dummy data
