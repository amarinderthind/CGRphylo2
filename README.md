# CGRphylo2: Chaos Game Representation for Phylogenetic Analysis

[![Bioconductor](https://img.shields.io/badge/Bioconductor-0.99.0-brightgreen.svg)](http://bioconductor.org/packages/CGRphylo2/)
[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview

**CGRphylo2** provides an efficient alignment-free approach for phylogenetic analysis of viral genomes using Chaos Game Representation (CGR). It accurately classifies closely related viral strains, including recombinants, making it particularly valuable during epidemic outbreaks where thousands of viral sequences need rapid analysis.

### About CGRphylo2

CGRphylo2 provides an alignment-free phylogenetic analysis method
    for viral genomes using Chaos Game Representation (CGR), a technique based
    on statistical physics concepts. Viruses exhibit high mutation rates,
    facilitating rapid evolution and emergence of new species, subspecies,
    strains, and recombinant forms. Accurate classification is crucial for
    understanding viral evolution and therapeutic development. Traditional
    phylogenetic methods require sequence alignment, which is computationally
    intensive. CGRphylo2 addresses this by implementing CGR-based whole-genome
    comparison that is fast, accurate, and computationally efficient. The
    package successfully classifies closely related viral lineages (demonstrated
    on SARS-CoV-2 lineages A and B), identifies recombinants (such as the XBB
    variant), and distinguishes multiple strains simultaneously. It processes
    sequences 5-13.7x faster than alignment-based methods (Clustal-Omega) with
    linear computational scaling. As a k-mer based approach, it enables
    simultaneous comparison of numerous closely-related sequences of different
    lengths. The package creates frequency matrices for distance calculations
    and phylogenetic tree construction, with outputs compatible with standard
    formats (MEGA, PHYLIP, Newick). Methods are based on Thind and Sinha (2023)
    <doi:10.2174/1389202924666230517115655>.
    
### Why CGRphylo?

- ✨ **Precision**: Accurately classifies closely related viral strains and recombinants
- ✨ **Speed**: 5-13.7× faster than traditional alignment methods (Clustal-Omega)
- ✨ **Scalability**: Linear computational cost with dataset size
- ✨ **Accessibility**: Designed for both high and low-resource settings

## Computational Efficiency

🚀 **CGRphylo2** processed 69 SARS-CoV-2 genomes **5 times faster** than Clustal-Omega.

🌐 For a dataset of 106 genomes, **CGRphylo2** outpaced Clustal-Omega by an incredible **13.7 times**.

In the world of MSAs (Multiple Sequence Alignment), computational costs skyrocket as datasets grow. Not for CGRphylo2! Adding a sequence is a breeze – just one frequency matrix calculation, breaking free from the computational intensity that others face.

## Installation

### From Bioconductor (Recommended)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CGRphylo2")
```

### Development Version

```r
# Install from GitHub
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("amarinderthind/CGRphylo2")
```

## Quick Start

```r
library(CGRphylo2)
library(seqinr)

# Load sequences
fastafile <- seqinr::read.fasta("sequences.fasta", 
                                seqtype = "DNA", 
                                as.string = TRUE)

# Filter sequences
fasta_filtered <- fastafile_new(fastafile, N_filter = 50)

# Calculate CGR frequency matrices (parallel)
freq_matrices <- parallelCGR(fasta_filtered, 
                             k_mer = 6, 
                             len_trim = min(sapply(fasta_filtered, nchar)))

# Calculate distance matrix
distance_matrix <- calculateDistanceMatrix(freq_matrices)

# Export results
saveMegaDistance("output.meg", distance_matrix)
savePhylipDistance("output.phy", distance_matrix, mode = "relaxed")

# Build tree
library(ape)
tree <- nj(as.dist(distance_matrix))
plot(tree)
```

## Key Features

### 1. Alignment-Free Approach

CGRphylo uses Chaos Game Representation to convert sequences into frequency matrices without requiring alignment:

- No position-by-position comparison needed
- Handles sequences of different lengths
- Robust to insertions/deletions

### 2. Efficient K-mer Analysis

- Flexible k-mer sizes (typically 4-8)
- Captures sequence composition patterns
- Scale-invariant representation

### 3. Multiple Distance Metrics

- Euclidean (default, recommended)
- Manhattan
- Squared Euclidean

### 4. Visualization Tools

```r
# Visualize CGR plot
cgr_coords <- cgrplot(1)
plot(cgr_coords[,1], cgr_coords[,2], 
     main = "CGR Plot", 
     cex = 0.2, pch = 4)
```

### 5. Export to Standard Formats

- MEGA format (.meg)
- PHYLIP format (.phy)
- Newick format (.nwk)
- Nexus format (.nex)

## Workflow

```
Input FASTA
    ↓
Filter Sequences
    ↓
Calculate CGR Frequencies (k-mers)
    ↓
Compute Distance Matrix
    ↓
Export / Build Tree
    ↓
Visualize Results
```

## Use Cases

### Viral Genomics

- SARS-CoV-2 variant analysis
- Recombinant detection
- Outbreak investigation
- Lineage classification

### Bacterial Genomics

- Strain typing
- Population structure
- Horizontal gene transfer detection

### Large-Scale Phylogenetics

- Rapid screening of thousands of genomes
- Real-time outbreak surveillance
- Resource-limited settings

## Performance

| Dataset Size | Clustal-Omega | CGRphylo | Speedup |
|--------------|---------------|----------|---------|
| 69 genomes   | 5 min         | 1 min    | 5×      |
| 106 genomes  | 13.7 min      | 1 min    | 13.7×   |

*Benchmarks on SARS-CoV-2 genomes (~30 kb each)*

## Documentation

- [Package Vignette](vignettes/CGRphylo_introduction.Rmd)
- [Function Reference](man/)
- [GitHub Repository](https://github.com/amarinderthind/CGRphylo2)

## Citation

If you use CGRphylo in your research, please cite:

> Thind AS, Sinha S (2023). Using Chaos-Game-Representation for Analysing the 
> SARS-CoV-2 Lineages, Newly Emerging Strains and Recombinants. *Current Genomics*, 
> 24(3). https://doi.org/10.2174/1389202924666230517115655

Full article: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10761335/

## Support

- [GitHub Issues](https://github.com/amarinderthind/CGRphylo2/issues)
- [Bioconductor Support Site](https://support.bioconductor.org/)

## License

GPL-3

## Author

**Amarinder Singh Thind**
- Email: amarinder.thind@gmail.com
- GitHub: [@amarinderthind](https://github.com/amarinderthind)

## Acknowledgments

We acknowledge the National Network for Mathematical and Computational Biology (NNMCB), DST, India for the internship programme at IISER Mohali for the initial part of the project.

## References

1. Thind AS, Sinha S (2023). Using Chaos-Game-Representation for Analysing the 
   SARS-CoV-2 Lineages. *Current Genomics*, 24(3).

2. Jeffrey HJ (1990). Chaos game representation of gene structure. 
   *Nucleic Acids Research*, 18(8):2163-2170.
