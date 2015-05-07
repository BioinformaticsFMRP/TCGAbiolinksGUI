library(stringr)
context("Search test")

test_that("u87 is mapped into nervous system", {
  suppressMessages(
  res <- biomics.search("u87")
  )
  expect_equal(solution, "BTO:0001484")
})

test_that("brain is mapped into nervous system", {
  suppressMessages(
    res <- biomics.search("brain")
  )
  expect_equal(solution, "BTO:0001484")
})

test_that("GBM is mapped into nervous system", {
  suppressMessages(
    res <- biomics.search("GBM")
  )
  expect_equal(solution, "BTO:0001484")
})

test_that("STAD is mapped into digestive system", {
  suppressMessages(
    res <- biomics.search("STAD")
  )
  expect_equal(solution, "BTO:0001491")
})

test_that("induced pluripotent stem cell is mapped into stem cell", {
  suppressMessages(
    res <- biomics.search("induced pluripotent stem cell")
  )
  expect_equal(solution, "BTO:0002666")
})

test_that("pancreas is mapped into disgestive and endocrine", {
  suppressMessages(
    res <- biomics.search("pancreas")
  )
  expect_equal(solution, "BTO:0001491,BTO:0001488")
})

