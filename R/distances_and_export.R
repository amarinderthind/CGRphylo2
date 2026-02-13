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
#' MEGA (Molecular Evolutionary Genetics Analysis) is widely used for phylogenetic
#' analysis. This function creates a distance matrix file compatible with MEGA's
#' format specifications.
#'
#' @examples
#' \dontrun{
#' # Create a distance matrix
#' dist_mat <- matrix(runif(9), nrow = 3)
#' rownames(dist_mat) <- colnames(dist_mat) <- c("Seq1", "Seq2", "Seq3")
#' 
#' # Save in MEGA format
#' saveMegaDistance("distances.meg", dist_mat)
#' }
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
    writeLines(paste("!Format DataType=Distance DataFormat=LowerLeft NTaxa=", 
                     n, ";", sep = ""), con)
    writeLines("", con)
    
    # Write sequence names
    writeLines("[", con)
    for (i in seq_len(n)) {
        writeLines(paste("#", seq_names[i], sep = ""), con)
    }
    writeLines("]", con)
    writeLines("", con)
    
    # Write distance matrix (lower triangular)
    for (i in seq_len(n)) {
        line <- paste(seq_names[i], "\t", sep = "")
        if (i > 1) {
            values <- sapply(1:(i-1), function(j) {
                format(distance_matrix[i, j], scientific = FALSE)
            })
            line <- paste(line, paste(values, collapse = "\t"), sep = "")
        }
        writeLines(line, con)
    }
    
    close(con)
    message(paste("Distance matrix saved to", filename))
}


#' Save distance matrix in PHYLIP format
#'
#' Exports a distance matrix to PHYLIP format for phylogenetic analysis with
#' various bioinformatics tools.
#'
#' @param filename Character. Output filename (typically .txt or .phy extension).
#' @param distance_matrix Numeric matrix. A square distance matrix with row and
#'   column names corresponding to sequence identifiers.
#' @param mode Character. PHYLIP format mode:
#'   \itemize{
#'     \item "original": Original PHYLIP format (10 character limit for names)
#'     \item "relaxed": Relaxed PHYLIP format (allows up to 250 characters)
#'   }
#'
#' @return NULL. Writes file as a side effect.
#'
#' @details
#' PHYLIP format is widely supported by phylogenetic software. The original format
#' limits sequence names to 10 characters, while the relaxed format allows longer
#' names. The relaxed format is recommended for modern applications.
#'
#' @references
#' PHYLIP format specification: http://www.phylo.org/index.php/help/relaxed_phylip
#'
#' @examples
#' \dontrun{
#' # Create a distance matrix
#' dist_mat <- matrix(runif(9), nrow = 3)
#' rownames(dist_mat) <- colnames(dist_mat) <- c("Seq1", "Seq2", "Seq3")
#' 
#' # Save in relaxed PHYLIP format
#' savePhylipDistance("distances.phy", dist_mat, mode = "relaxed")
#' }
#'
#' @export
savePhylipDistance <- function(filename, distance_matrix, mode = "relaxed") {
    n <- nrow(distance_matrix)
    seq_names <- rownames(distance_matrix)
    
    # Truncate names if original mode
    if (mode == "original") {
        seq_names <- substr(seq_names, 1, 10)
        seq_names <- sprintf("%-10s", seq_names)  # Left-align and pad to 10 chars
    }
    
    # Open file connection
    con <- file(filename, "w")
    
    # Write header (number of taxa)
    writeLines(paste("  ", n, sep = ""), con)
    
    # Write distance matrix
    for (i in seq_len(n)) {
        # Format distances with proper spacing
        distances <- sapply(seq_len(n), function(j) {
            format(round(distance_matrix[i, j], 6), width = 12, justify = "right")
        })
        
        line <- paste(seq_names[i], paste(distances, collapse = ""), sep = "  ")
        writeLines(line, con)
    }
    
    close(con)
    message(paste("Distance matrix saved to", filename, "in", mode, "PHYLIP format"))
}


#' Calculate pairwise distance matrix for multiple sequences
#'
#' Computes a full pairwise distance matrix from a list of CGR frequency matrices.
#'
#' @param freq_matrices List. A named list of frequency matrices, one per sequence,
#'   as returned by cgat().
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
#' \dontrun{
#' # Assuming freq_matrices is a list of frequency matrices
#' dist_matrix <- calculateDistanceMatrix(freq_matrices, 
#'                                        distance_type = "Euclidean")
#' }
#'
#' @export
calculateDistanceMatrix <- function(freq_matrices, distance_type = "Euclidean") {
    n <- length(freq_matrices)
    seq_names <- names(freq_matrices)
    
    # Initialize distance matrix
    distance <- matrix(0, nrow = n, ncol = n)
    rownames(distance) <- colnames(distance) <- seq_names
    
    # Calculate pairwise distances
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            distance[i, j] <- matrixDistance(freq_matrices[[i]], 
                                             freq_matrices[[j]], 
                                             distance_type = distance_type)
        }
    }
    
    # Round to reasonable precision
    distance <- round(distance, 6)
    
    return(distance)
}


#' Parallel computation of CGR frequency matrices
#'
#' Efficiently calculates CGR frequency matrices for multiple sequences using
#' parallel processing.
#'
#' @param sequences List. A list of DNA sequences (from fastafile_new or similar).
#' @param k_mer Integer. The k-mer size for frequency calculation.
#' @param len_trim Integer. Length to trim all sequences to.
#' @param num_cores Integer. Number of CPU cores to use. If NULL, uses
#'   detectCores() - 1.
#'
#' @return Named list. A list of frequency matrices, one per sequence.
#'
#' @details
#' This function uses parallel processing to speed up the calculation of CGR
#' frequency matrices for large datasets. It automatically detects available
#' cores and leaves one free to prevent system freezing.
#'
#' Note: On Windows, parallel processing falls back to sequential processing.
#'
#' @examples
#' \dontrun{
#' library(seqinr)
#' fasta <- read.fasta("sequences.fasta", seqtype = "DNA", as.string = TRUE)
#' filtered <- fastafile_new(fasta, N_filter = 50)
#' 
#' freq_mats <- parallelCGR(filtered, k_mer = 6, len_trim = 29000)
#' }
#'
#' @export
parallelCGR <- function(sequences, k_mer, len_trim, num_cores = NULL) {
    # Set up global environment for cgat function
    assign("fasta_filtered", sequences, envir = .GlobalEnv)
    
    # Determine number of cores
    if (is.null(num_cores)) {
        num_cores <- parallel::detectCores() - 1
        num_cores <- max(1, num_cores)  # Ensure at least 1 core
    }
    
    seq_names <- names(sequences)
    n_seq <- length(sequences)
    
    # Process function wrapper
    process_sequence <- function(n) {
        result <- cgat(k_mer, n, len_trim)
        return(result)
    }
    
    message(paste("Processing", n_seq, "sequences using", num_cores, "cores..."))
    
    # Use parallel processing on Unix-like systems, sequential on Windows
    if (.Platform$OS.type == "unix") {
        freq_matrices <- parallel::mclapply(seq_len(n_seq), 
                                           process_sequence, 
                                           mc.cores = num_cores)
    } else {
        # Windows fallback
        message("Windows detected: using sequential processing")
        freq_matrices <- lapply(seq_len(n_seq), process_sequence)
    }
    
    # Assign names
    names(freq_matrices) <- seq_names
    
    message("Processing complete!")
    
    return(freq_matrices)
}
