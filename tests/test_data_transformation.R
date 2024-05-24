#test data_transformation.R


library(testthat)
library(DESeq2)
library(limma)

source('../functions/data_transformation.R')
# Create some synthetic data to use in the tests
# Example count matrix
# Example count matrix with 100 genes and 50 samples
count_data <- matrix(sample(10:1000, 100 * 50, replace = TRUE), nrow = 100, ncol = 50)
colnames(count_data) <- paste0("Sample", 1:50)
rownames(count_data) <- paste0("Gene", 1:100)

# Example sample metadata with balanced conditions and batches
sample_data <- data.frame(
  X = colnames(count_data),
  Condition = factor(rep(c("Control", "Treatment"), each = 25)),
  Batch = factor(rep(c("Batch1", "Batch2"), each = 25)),
  tissue = rep(c("Tissue1", "Tissue2"), each = 25)  # For qsmooth
)
rownames(sample_data) <- colnames(count_data)

test_that("Test rlog transformation", {
  outfile <- tempfile()
  normalize_remove_vsd_rlog_batchef(count_data, sample_data, 'nobatch', outfile, 'rlog')
  
  # Read the output back in to check its contents
  result <- read.table(outfile, header = TRUE, sep = "\t", check.names = FALSE)
  expect_equal(dim(result), c(100, 50))
  expect_true(is.matrix(as.matrix(result)))
})

test_that("Test vsd transformation", {
  outfile <- tempfile()
  normalize_remove_vsd_rlog_batchef(count_data, sample_data, 'nobatch', outfile, 'vsd')
  
  result <- read.table(outfile, header = TRUE, sep = "\t", check.names = FALSE)
  expect_equal(dim(result), c(100, 50))
  expect_true(is.matrix(as.matrix(result)))
})

test_that("Test normalization", {
  outfile <- tempfile()
  normalize_remove_vsd_rlog_batchef(count_data, sample_data, 'nobatch', outfile, 'normalize')
  
  result <- read.table(outfile, header = TRUE, sep = "\t", check.names = FALSE)
  expect_equal(dim(result), c(100, 50))
  expect_true(is.matrix(as.matrix(result)))
})

test_that("Test quantilenorm_log2 transformation", {
  outfile <- tempfile()
  normalize_remove_vsd_rlog_batchef(count_data, sample_data, 'nobatch', outfile, 'quantilenorm_log2')
  
  result <- read.table(outfile, header = TRUE, sep = "\t", check.names = FALSE)
  expect_equal(dim(result), c(100, 50))
  expect_true(is.matrix(as.matrix(result)))
})

test_that("Test batch effect removal", {
  outfile <- tempfile()
  normalize_remove_vsd_rlog_batchef(count_data, sample_data, 'Batch', outfile, 'rlog')
  
  result <- read.table(outfile, header = TRUE, sep = "\t", check.names = FALSE)
  expect_equal(dim(result), c(100, 50))
  expect_true(is.matrix(as.matrix(result)))
})

