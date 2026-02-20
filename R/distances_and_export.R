#' Save distance matrix in MEGA format
#'
#' Exports a distance matrix to MEGA format for use with MEGA software for
#' phylogenetic tree visualization and analysis.
#'
#' @param filename Character. Output filename (typically with .meg extension).
#' @param distance_matrix Numeric matrix. A square distance matrix with row and
#'   column names corresponding to sequence identifiers.
#'
#' @return NULL. Writes file as a side effect.
#'
#' @details
#' MEGA (Molecular Evolutionary Genetics Analysis) is widely used for
#' phylogenetic analysis. This function creates a distance matrix file
#' compatible with MEGA's format specifications.
#'
#' @examples
#' # Build a small symmetric distance matrix
#' dist_mat <- matrix(
#'     c(0.00, 0.12, 0.25,
#'       0.12, 0.00, 0.18,
#'       0.25, 0.18, 0.00),
#'     nrow = 3
#' )
#' rownames(dist_mat) <- colnames(dist_mat) <- c("Seq1", "Seq2", "Seq3")
#'
#' # Save to a temporary file (no permanent files written during checks)
#' out <- tempfile(fileext = ".meg")
#' saveMegaDistance(out, dist_mat)
#'
#' # Inspect the first few lines of the output
#' writeLines(readLines(out, n = 8))
#'
#' @export
saveMegaDistance <- function(filename, distance_matrix) {
  n <- nrow(distance_matrix)
  seq_names <- rownames(distance_matrix)
  
  # Open file connection
  con <- file(filename, "w")
  
  # Write header
  writeLines("#mega", con)
  writeLines("!Title: Distance Matrix;", con)
  writeLines(paste0(
    "!Format DataType=Distance DataFormat=LowerLeft NTaxa=", n, ";"
  ), con)
  writeLines("", con)
  
  # Write sequence names
  writeLines("[", con)
  for (i in seq_len(n)) {
    writeLines(paste0("#", seq_names[i]), con)
  }
  writeLines("]", con)
  writeLines("", con)
  
  # Write distance matrix (lower triangular)
  for (i in seq_len(n)) {
    line <- paste0(seq_names[i], "\t")
    if (i > 1) {
      # FIX: was sapply(1:(i-1), ...); use vapply + seq_len
      values <- vapply(seq_len(i - 1), function(j) {
        format(distance_matrix[i, j], scientific = FALSE)
      }, FUN.VALUE = character(1))
      line <- paste0(line, paste(values, collapse = "\t"))
    }
    writeLines(line, con)
  }
  
  close(con)
  # FIX: was message(paste(...)); message() concatenates natively
  message("Distance matrix saved to ", filename)
}


#' Save distance matrix in PHYLIP format
#'
#' Exports a distance matrix to PHYLIP format for phylogenetic analysis with
#' various bioinformatics tools.
#'
#' @param filename Character. Output filename (typically .txt or .phy
#'   extension).
#' @param distance_matrix Numeric matrix. A square distance matrix with row
#'   and column names corresponding to sequence identifiers.
#' @param mode Character. PHYLIP format mode:
#'   \itemize{
#'     \item "original": Original PHYLIP format (10 character limit for names)
#'     \item "relaxed": Relaxed PHYLIP format (allows up to 250 characters)
#'   }
#'
#' @return NULL. Writes file as a side effect.
#'
#' @details
#' PHYLIP format is widely supported by phylogenetic software. The original
#' format limits sequence names to 10 characters, while the relaxed format
#' allows longer names. The relaxed format is recommended for modern
#' applications.
#'
#' @references
#' PHYLIP format specification:
#' http://www.phylo.org/index.php/help/relaxed_phylip
#'
#' @examples
#' # Build a small symmetric distance matrix
#' dist_mat <- matrix(
#'     c(0.00, 0.12, 0.25,
#'       0.12, 0.00, 0.18,
#'       0.25, 0.18, 0.00),
#'     nrow = 3
#' )
#' rownames(dist_mat) <- colnames(dist_mat) <- c("Seq1", "Seq2", "Seq3")
#'
#' # Save in relaxed PHYLIP format to a temporary file
#' out <- tempfile(fileext = ".phy")
#' savePhylipDistance(out, dist_mat, mode = "relaxed")
#'
#' # Inspect output
#' writeLines(readLines(out))
#'
#' @export
savePhylipDistance <- function(filename, distance_matrix, mode = "relaxed") {
  n <- nrow(distance_matrix)
  seq_names <- rownames(distance_matrix)
  
  # Truncate names if original mode
  if (mode == "original") {
    seq_names <- substr(seq_names, 1, 10)
    seq_names <- sprintf("%-10s", seq_names)
  }
  
  # Open file connection
  con <- file(filename, "w")
  
  # Write header (number of taxa)
  writeLines(paste0("  ", n), con)
  
  # Write distance matrix
  for (i in seq_len(n)) {
    # FIX: was sapply(); use vapply() to enforce return type
    distances <- vapply(seq_len(n), function(j) {
      format(round(distance_matrix[i, j], 6), width = 12,
             justify = "right")
    }, FUN.VALUE = character(1))
    
    line <- paste(seq_names[i], paste(distances, collapse = ""), sep = "  ")
    writeLines(line, con)
  }
  
  close(con)
  # FIX: was message(paste(...)); message() concatenates natively
  message("Distance matrix saved to ", filename, " in ", mode, " PHYLIP format")
}


#' Calculate pairwise distance matrix for multiple sequences
#'
#' Computes a full pairwise distance matrix from a list of CGR frequency
#' matrices.
#'
#' @param freq_matrices List. A named list of frequency matrices, one per
#'   sequence, as returned by cgat().
#' @param distance_type Character. Type of distance to calculate:
#'   "Euclidean" (default), "S_Euclidean", or "Manhattan".
#'
#' @return Numeric matrix. A symmetric distance matrix with sequence names as
#'   row and column names.
#'
#' @details
#' This function calculates all pairwise distances between sequences based on
#' their CGR frequency matrices. The resulting distance matrix can be used for
#' phylogenetic tree construction, clustering, or other downstream analyses.
#'
#' @examples
#' # Build two minimal frequency matrices (4 k-mers each, normalized)
#' fm1 <- matrix(c(0.4, 0.3, 0.2, 0.1), ncol = 1,
#'               dimnames = list(c("A", "C", "G", "T"), NULL))
#' fm2 <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 1,
#'               dimnames = list(c("A", "C", "G", "T"), NULL))
#' fm3 <- matrix(c(0.25, 0.25, 0.25, 0.25), ncol = 1,
#'               dimnames = list(c("A", "C", "G", "T"), NULL))
#'
#' freq_list <- list(Seq1 = fm1, Seq2 = fm2, Seq3 = fm3)
#'
#' dist_matrix <- calculateDistanceMatrix(freq_list, distance_type = "Euclidean")
#' print(round(dist_matrix, 4))
#'
#' @export
calculateDistanceMatrix <- function(freq_matrices,
                                    distance_type = "Euclidean") {
  n <- length(freq_matrices)
  seq_names <- names(freq_matrices)
  
  # Initialize distance matrix
  distance <- matrix(0, nrow = n, ncol = n)
  rownames(distance) <- colnames(distance) <- seq_names
  
  # Calculate pairwise distances
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      distance[i, j] <- matrixDistance(
        freq_matrices[[i]],
        freq_matrices[[j]],
        distance_type = distance_type
      )
    }
  }
  
  distance <- round(distance, 6)
  return(distance)
}


#' Parallel computation of CGR frequency matrices
#'
#' Efficiently calculates CGR frequency matrices for multiple sequences using
#' parallel processing.
#'
#' @param sequences List. A list of DNA sequences (from fastafile_new or
#'   similar).
#' @param k_mer Integer. The k-mer size for frequency calculation.
#' @param len_trim Integer. Length to trim all sequences to.
#' @param num_cores Integer. Number of CPU cores to use. If NULL, uses
#'   detectCores() - 1. Always set to 1 inside vignettes and examples.
#'
#' @return Named list. A list of frequency matrices, one per sequence.
#'
#' @details
#' This function uses parallel processing to speed up the calculation of CGR
#' frequency matrices for large datasets. It automatically detects available
#' cores and leaves one free to prevent system freezing.
#'
#' Note: On Windows, parallel processing falls back to sequential processing.
#' Always pass num_cores = 1 when calling from vignettes or examples to comply
#' with CRAN/Bioconductor check policies.
#'
#' @examples
#' # Small synthetic sequences — num_cores = 1 required in checked examples
#' seqs <- list(
#'     Seq1 = "ATCGATCGATCGATCGATCG",
#'     Seq2 = "GCTAGCTAGCTAGCTAGCTA"
#' )
#' assign("fasta_filtered", seqs, envir = .GlobalEnv)
#'
#' freq_mats <- parallelCGR(seqs, k_mer = 2, len_trim = 20, num_cores = 1)
#' cat("Matrices computed:", length(freq_mats), "\n")
#'
#' rm(fasta_filtered, envir = .GlobalEnv)
#'
#' @export
parallelCGR <- function(sequences, k_mer, len_trim, num_cores = NULL) {
  # Set up global environment for cgat function
  assign("fasta_filtered", sequences, envir = .GlobalEnv)
  
  # Determine number of cores
  if (is.null(num_cores)) {
    num_cores <- parallel::detectCores() - 1
    num_cores <- max(1, num_cores)
  }
  
  seq_names <- names(sequences)
  n_seq <- length(sequences)
  
  process_sequence <- function(n) {
    cgat(k_mer, n, len_trim)
  }
  
  # FIX: was message(paste(...)); message() concatenates natively
  message("Processing ", n_seq, " sequences using ", num_cores, " cores...")
  
  # Use parallel processing on Unix-like systems, sequential on Windows
  if (.Platform$OS.type == "unix") {
    freq_matrices <- parallel::mclapply(
      seq_len(n_seq),
      process_sequence,
      mc.cores = num_cores
    )
  } else {
    message("Windows detected: using sequential processing")
    freq_matrices <- lapply(seq_len(n_seq), process_sequence)
  }
  
  names(freq_matrices) <- seq_names
  message("Processing complete!")
  return(freq_matrices)
}