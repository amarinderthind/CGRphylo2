#' Calculate CGR frequency matrix for a DNA sequence
#'
#' This function computes the Chaos Game Representation (CGR) frequency matrix
#' for a given DNA sequence at a specified k-mer length. The CGR approach is an
#' alignment-free method that represents genomic sequences as frequency matrices,
#' which can then be used for phylogenetic analysis and sequence comparison.
#'
#' @param k_mer Integer. The word length (k-mer size) for frequency calculation.
#'   Typical values range from 4 to 8. Default is 6.
#' @param seq_index Integer. The index of the sequence in the global fasta_filtered
#'   object to process.
#' @param len_trim Integer. The length to which sequences should be trimmed for
#'   consistent comparison. Usually set to the minimum sequence length in the dataset.
#'
#' @return A matrix containing the frequencies of all possible k-mers for the
#'   given sequence. The matrix has 4^k rows corresponding to all possible k-mers.
#'
#' @details
#' The function creates a frequency matrix based on the CGR algorithm. Each k-mer's
#' frequency is calculated from the DNA sequence, and the resulting matrix can be
#' used to compute distances between sequences without requiring sequence alignment.
#'
#' The CGR method is particularly efficient for large datasets as the computational
#' cost of adding a new sequence is just one frequency matrix calculation, unlike
#' multiple sequence alignment where the cost increases quadratically.
#'
#' @examples
#' \dontrun{
#' # Assuming fasta_filtered is loaded
#' freq_matrix <- cgat(k_mer = 6, seq_index = 1, len_trim = 29000)
#' }
#'
#' @references
#' Thind AS, Sinha S (2023). Using Chaos-Game-Representation for Analysing the
#' SARS-CoV-2 Lineages, Newly Emerging Strains and Recombinants. Current Genomics,
#' 24(3). doi:10.2174/1389202924666230517115655
#'
#' @export
cgat <- function(k_mer, seq_index, len_trim) {
    # Get the sequence from the global fasta_filtered object
    sequence <- fasta_filtered[[seq_index]]
    
    # Trim sequence to specified length
    sequence <- substr(sequence, 1, len_trim)
    
    # Convert to uppercase
    sequence <- toupper(sequence)
    
    # Remove any non-ACGT characters
    sequence <- gsub("[^ACGT]", "", sequence)
    
    # Calculate total number of possible k-mers (4^k)
    n_kmers <- 4^k_mer
    
    # Generate all possible k-mers
    bases <- c("A", "C", "G", "T")
    all_kmers <- expand.grid(rep(list(bases), k_mer))
    all_kmers <- apply(all_kmers, 1, paste, collapse = "")
    
    # Initialize frequency matrix
    freq_matrix <- matrix(0, nrow = n_kmers, ncol = 1)
    rownames(freq_matrix) <- all_kmers
    
    # Count k-mer frequencies
    seq_length <- nchar(sequence)
    if (seq_length >= k_mer) {
        for (i in 1:(seq_length - k_mer + 1)) {
            kmer <- substr(sequence, i, i + k_mer - 1)
            if (kmer %in% all_kmers) {
                freq_matrix[kmer, 1] <- freq_matrix[kmer, 1] + 1
            }
        }
        
        # Normalize frequencies
        total_kmers <- sum(freq_matrix)
        if (total_kmers > 0) {
            freq_matrix <- freq_matrix / total_kmers
        }
    }
    
    return(freq_matrix)
}


#' Calculate distance between two CGR frequency matrices
#'
#' Computes the distance between two CGR frequency matrices using the specified
#' distance metric.
#'
#' @param matrix1 Numeric matrix. First CGR frequency matrix.
#' @param matrix2 Numeric matrix. Second CGR frequency matrix.
#' @param distance_type Character. Type of distance to calculate. Options are:
#'   \itemize{
#'     \item "Euclidean" (default): Standard Euclidean distance
#'     \item "S_Euclidean": Squared Euclidean distance
#'     \item "Manhattan": Manhattan (city block) distance
#'   }
#'
#' @return Numeric. The calculated distance between the two matrices.
#'
#' @details
#' This function calculates pairwise distances between CGR frequency matrices.
#' The Euclidean distance is most commonly used, but Manhattan and squared
#' Euclidean distances are also available for specific applications.
#'
#' @examples
#' \dontrun{
#' # Calculate Euclidean distance between two sequences
#' dist <- matrixDistance(freq_mat1, freq_mat2, distance_type = "Euclidean")
#' }
#'
#' @export
matrixDistance <- function(matrix1, matrix2, distance_type = "Euclidean") {
    if (distance_type == "Euclidean") {
        dist <- sqrt(sum((matrix1 - matrix2)^2))
    } else if (distance_type == "S_Euclidean") {
        dist <- sum((matrix1 - matrix2)^2)
    } else if (distance_type == "Manhattan") {
        dist <- sum(abs(matrix1 - matrix2))
    } else {
        stop("Invalid distance_type. Choose 'Euclidean', 'S_Euclidean', or 'Manhattan'")
    }
    return(dist)
}


#' Filter FASTA sequences by N content
#'
#' Filters a FASTA file by removing sequences with too many ambiguous (N) bases.
#'
#' @param fastafile List. A list of DNA sequences read by seqinr::read.fasta()
#' @param N_filter Integer. Maximum number of N bases allowed in a sequence.
#'   Sequences with more N's than this threshold will be removed.
#'
#' @return List. Filtered FASTA sequences.
#'
#' @details
#' This function is useful for quality control before phylogenetic analysis.
#' Sequences with excessive ambiguous bases can affect the accuracy of distance
#' calculations and tree construction.
#'
#' @examples
#' \dontrun{
#' library(seqinr)
#' fasta <- read.fasta("sequences.fasta", seqtype = "DNA", as.string = TRUE)
#' filtered_fasta <- fastafile_new(fasta, N_filter = 50)
#' }
#'
#' @export
fastafile_new <- function(fastafile, N_filter) {
    # Count N's in each sequence
    n_counts <- sapply(fastafile, function(seq) {
        length(grep("N", strsplit(toupper(seq), "")[[1]]))
    })
    
    # Filter sequences
    filtered <- fastafile[n_counts <= N_filter]
    
    message(paste("Filtered", length(fastafile) - length(filtered), 
                  "sequences with >", N_filter, "N bases"))
    
    return(filtered)
}


#' Create metadata table for sequences
#'
#' Extracts and compiles metadata including sequence length, GC content, and
#' N content from FASTA sequences.
#'
#' @param fastafile List. A list of DNA sequences read by seqinr::read.fasta()
#' @param N_filter Integer. N filter threshold (for reference in output)
#'
#' @return data.frame. A data frame with columns: name, length, GC_content, N_content
#'
#' @details
#' This function provides useful summary statistics for quality control and
#' understanding sequence characteristics before phylogenetic analysis.
#'
#' @examples
#' \dontrun{
#' library(seqinr)
#' fasta <- read.fasta("sequences.fasta", seqtype = "DNA", as.string = TRUE)
#' meta <- create_meta(fasta, N_filter = 50)
#' }
#'
#' @export
create_meta <- function(fastafile, N_filter) {
    meta_df <- data.frame(
        name = names(fastafile),
        length = sapply(fastafile, nchar),
        GC_content = sapply(fastafile, function(seq) {
            bases <- strsplit(toupper(seq), "")[[1]]
            gc_count <- sum(bases %in% c("G", "C"))
            gc_count / length(bases) * 100
        }),
        N_content = sapply(fastafile, function(seq) {
            bases <- strsplit(toupper(seq), "")[[1]]
            n_count <- sum(bases == "N")
            n_count
        }),
        stringsAsFactors = FALSE
    )
    
    return(meta_df)
}
