# Reproducibility

This directory contains `genome_info.csv` file that lists the RefSeq IDs of all the genomes we have used in our analysis.
This file contains the following columns:

- `Genome`: The RefSeq ID of the genome
- `Kingdom`: The kingdom of the genome (e.g. Bacteria, Archaea)
- `Genome Length`: The length of the genome in base pairs
- `Unique k-mers`: The number of unique $k$-mers in the genome
- `k-mer Size`: The length of the $k$-mers used in the analysis
- `Alpha`: $\alpha^\ast$ value of `CaPLa` for the genome
- `Beta Min`: $\beta_\mathrm{low}^\ast$ value of `CaPLa` for the genome
- `Beta Max`: $\beta_\mathrm{high}^\ast$ value of `CaPLa` for the genome
- `Eps_{1, 2, ..., 1024}`: One column for each of the error bounds $\varepsilon \in \lbrace 1, 2, \ldots, 1024 \rbrace$ showing the number of segments in the PLA
