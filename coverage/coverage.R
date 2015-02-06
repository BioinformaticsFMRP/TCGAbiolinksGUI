loc <- file.path("/Users/tiago/prime")

require(testthat)

# example saturate - testthat partial coverage
reportCoverage(sourcefiles = list.files(file.path(loc, "R"), full.names = TRUE),
               executionfiles = list.files(file.path(loc, 
                                                      "tests", "testthat"), full.names = TRUE), 
               reportfile = "testCoverage_saturate_prime.html", 
               writereport = TRUE, clean = TRUE)