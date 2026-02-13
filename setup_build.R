#!/usr/bin/env Rscript
# CGRphylo Package Setup and Build Script
# Run this script to set up, document, build, and check the package

cat("==================================================\n")
cat("CGRphylo Package Setup and Build\n")
cat("==================================================\n\n")

# Install required packages
cat("Step 1: Installing required packages...\n")
required_packages <- c("devtools", "roxygen2", "testthat", "knitr", 
                      "rmarkdown", "BiocCheck")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    if (pkg == "BiocCheck") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("BiocCheck")
    } else {
      install.packages(pkg)
    }
  }
}

# Set working directory to package root
if (basename(getwd()) != "CGRphylo_bioconductor") {
  if (dir.exists("CGRphylo_bioconductor")) {
    setwd("CGRphylo_bioconductor")
  }
}

cat("\nStep 2: Generating documentation with roxygen2...\n")
devtools::document()

cat("\nStep 3: Building package...\n")
pkg_file <- devtools::build()

cat("\nStep 4: Running R CMD check...\n")
check_results <- devtools::check()

cat("\n==================================================\n")
cat("Check Results Summary:\n")
cat("==================================================\n")
cat(paste("Errors:", length(check_results$errors), "\n"))
cat(paste("Warnings:", length(check_results$warnings), "\n"))
cat(paste("Notes:", length(check_results$notes), "\n"))

if (length(check_results$errors) > 0) {
  cat("\nErrors found:\n")
  print(check_results$errors)
}

if (length(check_results$warnings) > 0) {
  cat("\nWarnings found:\n")
  print(check_results$warnings)
}

cat("\nStep 5: Running BiocCheck...\n")
BiocCheck::BiocCheck(".")

cat("\n==================================================\n")
cat("Setup Complete!\n")
cat("==================================================\n")
cat("\nNext steps:\n")
cat("1. Review any errors, warnings, or notes above\n")
cat("2. Run tests: devtools::test()\n")
cat("3. Build vignette: devtools::build_vignettes()\n")
cat("4. Install locally: devtools::install()\n")
cat("5. When ready, follow BIOCONDUCTOR_SUBMISSION.md\n")
cat("\n")
