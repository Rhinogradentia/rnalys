
library(testthat)
# Load your script
source("../functions/run_edgeR.R")

test_that("run_DE handles input without batch effects correctly", {
  gene_counts <- matrix(rnbinom(1000, mu=100, size=1), nrow=100, ncol=50)
  sample_info <- data.frame(condition = factor(rep(c("Control", "Treatment"), each=25)))
  rownames(sample_info) <- paste0("Sample", 1:50)
  
  design_formula <- "~condition"
  temp_output_file <- tempfile()
  
  run_DE(indata = gene_counts, insample = sample_info, design = design_formula,
         reference = "Control", outfile = temp_output_file)
  
  expect_true(file.exists(temp_output_file), info = "Output file should exist.")
  expect_true(file.info(temp_output_file)$size > 0, info = "Output file should not be empty.")
})

test_that("run_DE handles input with batch effects correctly", {
  gene_counts <- matrix(rnbinom(1000, mu=120, size=1), nrow=100, ncol=10)
  sample_info <- data.frame(Treatment = factor(rep(c("Control", "Treatment"), each=5)),
                            Batch = factor(rep(c("Batch1", "Batch2"), times=5)))
  rownames(sample_info) <- paste0("Sample", 1:10)
  
  design_formula <- "~ Batch + Treatment"
  temp_output_file <- tempfile()
  
  run_DE(indata = gene_counts, insample = sample_info, design = design_formula,
         reference = "Control", outfile = temp_output_file)
  
  expect_true(file.exists(temp_output_file), info = "Output file should exist.")
  expect_true(file.info(temp_output_file)$size > 0, info = "Output file should not be empty.")
})
