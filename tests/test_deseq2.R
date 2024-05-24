
library(testthat)
# Load your script
#rint(getwd())
source("../functions/run_deseq2_working.R")

# Testing extract_first_variable
test_that("extract_last_variable correctly extracts the last variable", {
  expect_equal(extract_last_variable("~condition"), "condition")
  expect_equal(extract_last_variable("~condition + treatment"), "treatment")
  expect_error(extract_last_variable("~"))
})

# Testing run_DE function with simulated data
test_that("run_DE handles basic input without error", {
  # Simulate gene counts data
  set.seed(123)
  gene_counts <- matrix(rnbinom(1000, mu=100, size=1), nrow=100, ncol=10)
  
  # Create sample information
  sample_info <- data.frame(
    condition = factor(rep(c("Control", "Treatment"), each=5))
  )
  rownames(sample_info) <- paste0("Sample", 1:10)

  # Compute row means (not directly used in this context)
  rowm <- 10
  
  # Define the design formula
  design_formula <- "~condition"
  
  # Specify a temporary file for output
  temp_output_file <- tempfile()

  #expect_silent(run_DE(gene_counts, sample_info, rowm, design_formula, temp_output_file, reference="Control"))
  run_DE(gene_counts, sample_info, rowm, design_formula, temp_output_file, reference="Control")
  expect_true(file.exists(temp_output_file), info = "Output file should exist.")
  expect_true(file.info(temp_output_file)$size > 0, info = "Output file should not be empty.")

  #Checking results
  results <- as.data.frame(read.csv(temp_output_file, sep="\t"))
  expect_true(nrow(results) > 0, info = "Results should have non-zero rows")
  

})

#test when batch is present

test_that("run_DE handles batch effects without error", {
  # Simulate gene counts data
  set.seed(124)
  gene_counts <- matrix(rnbinom(1000, mu = 120, size = 1), nrow = 100, ncol = 50)
  
  # Create sample information including batch effect
  sample_info <- data.frame(
    Treatment = factor(rep(c("Control", "Treatment"), each = 25)),
    Batch = factor(rep(c("Batch1", "Batch2"), times = 25))
  )
  rownames(sample_info) <- paste0("Sample", 1:50)
  # Updated design formula to include batch effect
  design_formula <- "~ Batch + Treatment"
  
  # Temporary file for output
  temp_output_file <- tempfile()
  
  # Run the DE analysis
  run_DE(
    indata = gene_counts,
    insample = sample_info,
    rowm = rowMeans(gene_counts),
    design = design_formula,
    outfile = temp_output_file,
    reference = "Control"
  )
  
  # Check if the output file was created and is not empty
  expect_true(file.exists(temp_output_file), info = "Output file should exist.")
  expect_true(file.info(temp_output_file)$size > 0, info = "Output file should not be empty.")
  
  # Optionally read and check results content
  results <- read.csv(temp_output_file)
  #print(results)
  
  # Ensure results dataframe is correctly structured
  expect_true(nrow(results) > 0, info = "Results should have non-zero rows")
  #expect_true(ncol(results) >= 5, info = "Results should have at least 5 columns for DESeq2 output")
})

# Run all tests
