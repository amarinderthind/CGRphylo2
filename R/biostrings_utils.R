#' Convert a DNAStringSet to a named list of sequences
#'
#' Converts a \code{\link[Biostrings]{DNAStringSet}} object to a named list of
#' character strings suitable for use with CGRphylo2 functions such as
#' \code{fastafile_new}, \code{create_meta}, and \code{parallelCGR}.
#'
#' @param dna A \code{DNAStringSet} object containing one or more DNA sequences.
#'
#' @return A named list of character strings, one element per sequence.
#'
#' @details
#' Bioconductor workflows commonly store DNA sequences as
#' \code{DNAStringSet} objects. This function bridges that format with
#' CGRphylo2's internal list representation, allowing seamless use of
#' Bioconductor data structures in CGR-based phylogenetic analysis.
#'
#' @examples
#' if (requireNamespace("Biostrings", quietly = TRUE)) {
#'     dna <- Biostrings::DNAStringSet(c(
#'         seq1 = "ATCGATCGATCGATCG",
#'         seq2 = "GCTAGCTAGCTAGCTA"
#'     ))
#'     seqs <- from_DNAStringSet(dna)
#'     length(seqs)   # 2
#'     nchar(seqs[[1]])  # 16
#' }
#'
#' @importFrom Biostrings DNAStringSet
#' @export
from_DNAStringSet <- function(dna) {
    seqs <- as.list(as.character(dna))
    names(seqs) <- names(dna)
    seqs
}
