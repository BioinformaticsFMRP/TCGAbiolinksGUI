# Load a unit testing framework and testCoverage.
require(testthat)
require("testCoverage")
loc <- file.path(".")
# example saturate - testthat partial coverage
reportCoverage(sourcefiles = list.files(file.path(loc, "R"),
                                        full.names = TRUE, pattern = ".R"),
               executionfiles = list.files(
                    file.path(loc,"tests", "testthat"),
                    full.names = TRUE
               ),
               reportfile = "R/coverage/test_coverage.html",
               writereport = TRUE, clean = TRUE)
