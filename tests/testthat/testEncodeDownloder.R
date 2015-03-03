library(testthat)
context("Verify Parameters")

# TARGET
test_that("Parameter target should have the prefix target.investigated_as=  ", {
  expect_equal(formatTarget("transcription factor"),
               "target.investigated_as=transcription%20factor"
               )
})

test_that("If more than one parameter target,
          both should have the prefix target.investigated_as=  ", {
  expect_equal(formatTarget(c("transcription factor","tag")),
               paste ("target.investigated_as=transcription%20factor&",
                      "target.investigated_as=tag",
                      sep = ""
               ))
})

test_that("A null parameter returns NULL", {
  expect_equal(formatTarget(NULL),NULL)
})

# Sample
test_that("Parameter sample should have the
          prefix replicates.library.biosample.biosample_type=  ", {
  expect_equal(formatSample("primary cell"),
               "replicates.library.biosample.biosample_type=primary%20cell"
               )
})

test_that("If more than one parameter sample, both should have the prefix  ", {
  expect_equal(formatSample(c("primary cell","tissue")),
               paste ("replicates.library.biosample.",
                      "biosample_type=primary%20cell&",
                      "replicates.library.biosample.",
                      "biosample_type=tissue",
                      sep = ""
               ))
})

test_that("A null parameter returns NULL", {
  expect_equal(formatSample(NULL),NULL)
})

# FType
test_that("Parameter fType should have the prefix files.file_format=  ", {
  expect_equal(formatFType("bam"), "files.file_format=bam")
})

test_that("If more than one parameter fType, both should have the prefix", {
  expect_equal(formatFType(c("bam","bigWig")),
               "files.file_format=bam&files.file_format=bigWig")
})

test_that("A null parameter returns NULL", {
  expect_equal(formatFType(NULL),NULL)
})

# Assay
test_that("Parameter assay should have the prefix assay_term_name=  ", {
  expect_equal(formatAssay("ChIP-seq"), "assay_term_name=ChIP-seq")
})

test_that("If more than one parameter assay,
           both should have the prefix assay_term_name=  ", {
  expect_equal(formatAssay( c("ChIP-seq", "RIP-chip")),
                "assay_term_name=ChIP-seq&assay_term_name=RIP-chip"
               )
})

test_that("A null parameter returns NULL", {
  expect_equal(formatAssay(NULL),NULL)
})


# Assembly
test_that("Parameter assembly should have the prefix assembly=  ", {
  expect_equal(formatAssembly("mm9"), "assembly=mm9")
})

test_that("If more than one parameter assembly,
          both should have the prefix assembly=  ", {
  expect_equal(formatAssembly(c("mm9","hg19")),
               "assembly=mm9&assembly=hg19"
               )
})

test_that("A null parameter returns NULL", {
  expect_equal(formatAssembly(NULL),NULL)
})


# Search
test_that("Parameter assembly should have the prefix searchTerm=  ", {
  expect_equal(formatSearch("bone"), "searchTerm=bone")
})

test_that("If more than one word in search, both should be separated by +  ", {
  expect_equal(formatSearch("bone chip"),
             "searchTerm=bone+chip")
})

test_that("A null parameter returns NULL", {
  expect_equal(formatSearch(NULL),NULL)
})

# Type
test_that("Parameter type should have the prefix type=  ", {
  expect_equal(formatType("experiment"), "type=experiment")
})


test_that("A null parameter returns NULL", {
  expect_equal(formatType(NULL),NULL)
})
