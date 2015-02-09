context("Return verification")

test_that("Prime of a prime number returns TRUE", {
  expect_true(prime(3))
  expect_true(prime(13))
  expect_true(prime(2))
})

test_that("Prime of a non-prime number returns FALSE", {
  expect_false(prime(220))
  expect_false(prime(9))
})

context("Handling NA")
test_that("Prime of missing is missing", {
  expect_error(prime(NA),cat("Error in prime(NA) : x must be a number","\n"))
})
