# Unit tests for CGRphylo2 package

library(testthat)
library(CGRphylo2)

# Test data
test_sequences <- list(
    seq1 = "ATCGATCGATCG",
    seq2 = "GCTAGCTAGCTA",
    seq3 = "AAAATTTTCCCCGGGG"
)

test_that("fastafile_new filters sequences correctly", {
    # Create test data with N's
    test_seqs_with_n <- list(
        good_seq = "ATCGATCG",
        bad_seq = "ATCGNNNNNATCG",  # 5 N's
        okay_seq = "ATCGNNATCG"      # 2 N's
    )
    
    # Filter with threshold of 3 N's
    filtered <- fastafile_new(test_seqs_with_n, N_filter = 3)
    
    expect_equal(length(filtered), 2)
    expect_true("good_seq" %in% names(filtered))
    expect_true("okay_seq" %in% names(filtered))
    expect_false("bad_seq" %in% names(filtered))
})

test_that("create_meta generates correct metadata", {
    meta <- create_meta(test_sequences, N_filter = 0)
    
    expect_s3_class(meta, "data.frame")
    expect_equal(nrow(meta), 3)
    expect_true(all(c("name", "length", "GC_content", "N_content") %in% colnames(meta)))
    
    # Check sequence lengths
    expect_equal(meta$length[1], nchar(test_sequences$seq1))
    expect_equal(meta$length[2], nchar(test_sequences$seq2))
    
    # Check GC content for seq3 (50% GC)
    expect_equal(meta$GC_content[3], 50)
})

test_that("cgat generates frequency matrix with correct dimensions", {
    # Set up global variable
    assign("fasta_filtered", test_sequences, envir = .GlobalEnv)
    
    k_mer <- 3
    len_trim <- 12
    
    freq_matrix <- cgat(k_mer = k_mer, seq_index = 1, len_trim = len_trim)
    
    # Check dimensions (4^3 = 64 possible 3-mers)
    expect_equal(nrow(freq_matrix), 4^k_mer)
    expect_equal(ncol(freq_matrix), 1)
    
    # Check that frequencies sum to 1 (normalized)
    expect_equal(sum(freq_matrix), 1, tolerance = 1e-6)
    
    # All frequencies should be non-negative
    expect_true(all(freq_matrix >= 0))
})

test_that("matrixDistance calculates distances correctly", {
    # Create two identical matrices
    mat1 <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 1)
    mat2 <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 1)
    mat3 <- matrix(c(0.2, 0.3, 0.4, 0.5), ncol = 1)
    
    # Distance between identical matrices should be 0
    dist_euclidean <- matrixDistance(mat1, mat2, distance_type = "Euclidean")
    expect_equal(dist_euclidean, 0)
    
    # Test Manhattan distance
    dist_manhattan <- matrixDistance(mat1, mat3, distance_type = "Manhattan")
    expect_equal(dist_manhattan, 0.4)  # |0.1-0.2| + |0.2-0.3| + |0.3-0.4| + |0.4-0.5|
    
    # Test Squared Euclidean
    dist_sq_euclidean <- matrixDistance(mat1, mat3, distance_type = "S_Euclidean")
    expect_equal(dist_sq_euclidean, 0.04)  # 4 * (0.1)^2
})

test_that("cgrplot generates coordinates correctly", {
    assign("fasta_filtered", test_sequences, envir = .GlobalEnv)
    
    cgr_coords <- cgrplot(seq_index = 1)
    
    # Check output is a matrix
    expect_true(is.matrix(cgr_coords))
    expect_equal(ncol(cgr_coords), 2)
    
    # Check column names
    expect_equal(colnames(cgr_coords), c("x", "y"))
    
    # All coordinates should be between 0 and 1
    expect_true(all(cgr_coords >= 0 & cgr_coords <= 1))
})

test_that("calculateDistanceMatrix produces symmetric matrix", {
    assign("fasta_filtered", test_sequences, envir = .GlobalEnv)
    
    # Create frequency matrices
    freq_matrices <- list(
        seq1 = cgat(3, 1, 12),
        seq2 = cgat(3, 2, 12),
        seq3 = cgat(3, 3, 16)
    )
    
    dist_matrix <- calculateDistanceMatrix(freq_matrices, 
                                          distance_type = "Euclidean")
    
    # Check dimensions
    expect_equal(nrow(dist_matrix), 3)
    expect_equal(ncol(dist_matrix), 3)
    
    # Check symmetry
    expect_equal(dist_matrix, t(dist_matrix))
    
    # Diagonal should be zero
    expect_equal(diag(dist_matrix), c(0, 0, 0))
    
    # All values should be non-negative
    expect_true(all(dist_matrix >= 0))
})

test_that("saveMegaDistance creates valid file", {
    # Create a small test matrix
    test_matrix <- matrix(c(0, 0.1, 0.2,
                           0.1, 0, 0.15,
                           0.2, 0.15, 0), nrow = 3)
    rownames(test_matrix) <- colnames(test_matrix) <- c("A", "B", "C")
    
    temp_file <- tempfile(fileext = ".meg")
    saveMegaDistance(temp_file, test_matrix)
    
    # Check file exists
    expect_true(file.exists(temp_file))
    
    # Read and check content
    content <- readLines(temp_file)
    expect_true(any(grepl("#mega", content)))
    expect_true(any(grepl("NTaxa=3", content)))
    
    # Clean up
    unlink(temp_file)
})

test_that("savePhylipDistance creates valid file", {
    # Create a small test matrix
    test_matrix <- matrix(c(0, 0.1, 0.2,
                           0.1, 0, 0.15,
                           0.2, 0.15, 0), nrow = 3)
    rownames(test_matrix) <- colnames(test_matrix) <- c("Seq1", "Seq2", "Seq3")
    
    temp_file <- tempfile(fileext = ".phy")
    savePhylipDistance(temp_file, test_matrix, mode = "relaxed")
    
    # Check file exists
    expect_true(file.exists(temp_file))
    
    # Read and check content
    content <- readLines(temp_file)
    expect_true(grepl("^\\s*3", content[1]))  # First line should have number of taxa
    
    # Clean up
    unlink(temp_file)
})

test_that("Invalid distance types throw errors", {
    mat1 <- matrix(c(0.1, 0.2), ncol = 1)
    mat2 <- matrix(c(0.3, 0.4), ncol = 1)
    
    expect_error(matrixDistance(mat1, mat2, distance_type = "InvalidType"))
})

test_that("parallelCGR works correctly", {
    skip_on_cran()  # Skip on CRAN as parallel processing can be unstable
    
    freq_matrices <- parallelCGR(test_sequences, 
                                 k_mer = 3, 
                                 len_trim = 12,
                                 num_cores = 1)
    
    expect_equal(length(freq_matrices), 3)
    expect_equal(names(freq_matrices), names(test_sequences))
    
    # Each frequency matrix should have correct dimensions
    lapply(freq_matrices, function(mat) {
        expect_equal(nrow(mat), 4^3)
    })
})
