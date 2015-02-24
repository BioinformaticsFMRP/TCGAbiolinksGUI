library(testthat)

context("Roadmap download")

test_that("A compressed file is decompressed", {
  filename <- "roadmap_test.gz"
  df1 <- data.frame(id = seq(1,10,1), var1 = runif(10), var2 = runif(10))
  gz1 <- gzfile(filename, "w")
  write.csv(df1, gz1)
  close(gz1)
  expect_equal(.uncompress(filename)[[1]][1], "roadmap_test")
  file.remove("roadmap_test")
  file.remove("roadmap_test.gz")
})

test_that("If list is empty, return NULL", {
  expect_equal(.uncompress(NULL),NULL)
})

test_that("If no file is found, return empty list", {
  filename <- "roadmap_test2.gz"
  expect_equal(.uncompress(filename),list(NULL))
})

test_that("Test list of files", {
  filename <- "roadmap_test.gz"
  df1 <- data.frame(id = seq(1,10,1), var1 = runif(10), var2 = runif(10))
  gz1 <- gzfile(filename, "w")
  write.csv(df1, gz1)
  close(gz1)
  ret <- .uncompress( c("roadmap_test2.gz","roadmap_test.gz"))
  expect_equal(ret[[1]], NULL)
  expect_equal(ret[[2]][1], "roadmap_test")
  file.remove("roadmap_test")
  file.remove("roadmap_test.gz")
})


test_that("A directory is created if it does not exists", {
  if (file.exists("roadmap_test"))
      file.remove("roadmap_test")
  expect_true(.mkdir("roadmap_test"))
})

test_that("A directory is not created if it exists", {
  expect_equal(.mkdir("roadmap_test"), NULL)
  if (file.exists("roadmap_test"))
      file.remove("roadmap_test")
})
