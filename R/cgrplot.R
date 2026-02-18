# Avoid NOTE about global variables
utils::globalVariables(c("fasta_filtered"))

#' Generate CGR plot coordinates for a DNA sequence
#'
#' Creates x and y coordinates for visualizing a Chaos Game Representation (CGR)
#' plot of a DNA sequence. The CGR plot is a 2D representation that reveals
#' patterns and composition of genomic sequences.
#'
#' @param seq_index Integer. The index of the sequence in the global
#'   fasta_filtered object to plot.
#'
#' @return Matrix with two columns (x and y coordinates) for plotting the CGR.
#'   Each row represents the position of one nucleotide in the CGR space.
#'
#' @details
#' The Chaos Game Representation is an iterative mapping technique that creates
#' a 2D representation of DNA sequences. Each base (A, C, G, T) is assigned to
#' a corner of a unit square, and the sequence is plotted by iteratively moving
#' halfway from the current position to the corner corresponding to the next base.
#'
#' The resulting plot has fractal properties and reveals sequence composition,
#' repeats, and other genomic features. Different sequences create distinct
#' patterns that reflect their underlying genomic structure.
#'
#' Corner assignments (standard):
#' \itemize{
#'   \item A: (0, 0) - bottom left
#'   \item C: (1, 0) - bottom right
#'   \item G: (1, 1) - top right
#'   \item T: (0, 1) - top left
#' }
#'
#' @examples
#' \dontrun{
#' # Generate CGR coordinates for first sequence
#' cgr_coords <- cgrplot(1)
#'
#' # Plot the CGR
#' plot(cgr_coords[, 1], cgr_coords[, 2],
#'   main = "CGR Plot",
#'   xlab = "", ylab = "",
#'   cex = 0.2, pch = 4
#' )
#' }
#'
#' @references
#' Jeffrey HJ (1990). Chaos game representation of gene structure.
#' Nucleic Acids Research, 18(8):2163-2170.
#'
#' @export
cgrplot <- function(seq_index) {
  # Get the sequence from the global fasta_filtered object
  sequence <- fasta_filtered[[seq_index]]

  # Convert to uppercase and split into individual bases
  sequence <- toupper(sequence)
  bases <- strsplit(sequence, "")[[1]]

  # Remove non-ACGT bases
  bases <- bases[bases %in% c("A", "C", "G", "T")]

  # Define corner positions for each base
  corners <- list(
    A = c(0, 0),
    C = c(1, 0),
    G = c(1, 1),
    T = c(0, 1)
  )

  # Initialize starting position (center)
  current_pos <- c(0.5, 0.5)

  # Initialize matrix to store coordinates
  n_bases <- length(bases)
  cgr_coords <- matrix(0, nrow = n_bases, ncol = 2)
  colnames(cgr_coords) <- c("x", "y")

  # Calculate CGR coordinates
  for (i in seq_len(n_bases)) {
    base <- bases[i]
    if (base %in% names(corners)) {
      # Move halfway to the corner of the current base
      corner <- corners[[base]]
      current_pos <- (current_pos + corner) / 2
      cgr_coords[i, ] <- current_pos
    }
  }

  return(cgr_coords)
}


#' Plot CGR with customization options
#'
#' A convenience wrapper function to create a CGR plot with sensible defaults.
#'
#' @param seq_index Integer. The index of the sequence to plot.
#' @param main Character. Main title for the plot. If NULL, uses sequence name.
#' @param cex Numeric. Point size (default 0.2 for dense sequences).
#' @param pch Integer. Point character type (default 4 for crosses).
#' @param col Character. Color for points (default "black").
#' @param ... Additional arguments passed to plot().
#'
#' @return NULL. Creates a plot as a side effect.
#'
#' @examples
#' \dontrun{
#' # Simple CGR plot
#' plot_cgr(1)
#'
#' # Customized plot
#' plot_cgr(1, main = "My Sequence", col = "blue", cex = 0.3)
#' }
#'
#' @export
plot_cgr <- function(seq_index, main = NULL, cex = 0.2, pch = 4,
                     col = "black", ...) {
  cgr_coords <- cgrplot(seq_index)

  if (is.null(main)) {
    main <- paste("CGR plot of", names(fasta_filtered)[seq_index])
  }

  plot(cgr_coords[, 1], cgr_coords[, 2],
    main = main,
    xlab = "", ylab = "",
    cex = cex, pch = pch,
    col = col,
    frame.plot = TRUE,
    ...
  )
}
