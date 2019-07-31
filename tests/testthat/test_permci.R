context("Testing permci")

test_that("permci parallel works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)
  tmp1 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                     nperm = 10, ncores = 1, quietly = T)
  expect_true(getDoParWorkers() == 1)
  tmp2 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                     nperm = 10, ncores = 2, quietly = T)
  expect_true(getDoParWorkers() == 2)
  tmp3 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                     nperm = 10, ncores = 1, quietly = T)
  expect_true(getDoParWorkers() == 1)
})

test_that("permci seed parameter works as it should", {
  ds <- gendata_crt(nclus = c(5, 5), size = c(10, 10), theta = 0, sigma = 1,
                    mu = 0, sd = 1)

  # not parallel
  tmp1 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = 1, seed = 123, quietly = T)
  tmp2 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = 1, seed = 123, quietly = T)
  tmp3 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = 1, quietly = T)
  tmp4 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = 1, seed = 124, quietly = T)
  expect_true(identical(tmp1, tmp2))
  expect_false(identical(tmp2, tmp3))
  expect_false(identical(tmp2, tmp4))

  # parallel
  ncores <- parallel::detectCores() - 1
  tmp1 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = ncores, seed = 123, quietly = T)
  tmp2 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = ncores, seed = 123, quietly = T)
  tmp3 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = ncores, quietly = T)
  tmp4 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                       nperm = 100, ncores = ncores, seed = 124, quietly = T)
  expect_true(identical(tmp1, tmp2))
  expect_false(identical(tmp2, tmp3))
  expect_false(identical(tmp2, tmp4))

  # not parallel vs parallel
  tmp1 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                     nperm = 100, ncores = 1, seed = 123, quietly = T)
  tmp2 <- permci_glm(y ~ trt, trtname = 'trt', runit = 'clusid', data = ds,
                     nperm = 100, ncores = 2, seed = 123, quietly = T)
  expect_true(identical(tmp1$trace, tmp2$trace))
})

# should add other permci_xxx tests, but first need dummy data
