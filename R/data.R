#' Example SARS-CoV-2 sequences for testing
#'
#' A small dataset of SARS-CoV-2 genome sequences for demonstrating
#' CGRphylo functionality.
#'
#' @format A list with 5 DNA sequences:
#' \describe{
#'   \item{sequence_1}{First SARS-CoV-2 genome sequence}
#'   \item{sequence_2}{Second SARS-CoV-2 genome sequence}
#'   \item{sequence_3}{Third SARS-CoV-2 genome sequence}
#'   \item{sequence_4}{Fourth SARS-CoV-2 genome sequence}
#'   \item{sequence_5}{Fifth SARS-CoV-2 genome sequence}
#' }
#'
#' @details
#' These sequences are truncated versions (~5kb) of actual SARS-CoV-2 genomes
#' for demonstration purposes. They represent different lineages and can be
#' used to test the CGRphylo pipeline.
#'
#' @source Simulated data based on NCBI SARS-CoV-2 sequences
#'
#' @examples
#' # This is a placeholder - actual data would be in inst/extdata
#' # To use real example data:
#' # example_file <- system.file("extdata", "example_seqs.fasta", 
#' #                             package = "CGRphylo2")
#' # sequences <- seqinr::read.fasta(example_file, seqtype = "DNA")
#'
"example_sequences"
