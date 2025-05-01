# CaPLa

The **Canonical Piecewise Linear approximability (CaPLa)** is a measure for quantifying the efficiency of learned data structures based on piecewise linear approximations (PLAs) on a given dataset.

CaPLa captures the *tightest power-law bound* on the average number of elements spanned by each segment in the PLA.
It is defined as a triple  $(\alpha^\ast, \beta_\mathrm{low}^\ast, \beta_\mathrm{high}^\ast)$ such that,
for a dataset with $n$ distinct elements and a PLA using $b(\varepsilon)$ segments at an error bound $\varepsilon$, the following inequality holds for all $\varepsilon$ in a given set of error bounds:
$$\beta_\mathrm{low}^\ast \varepsilon^\alpha \leq \frac{n}{b(\varepsilon)} \leq \beta_\mathrm{high}^\ast \varepsilon^\alpha.$$

This repository provides a **tool for computing the CaPLa of a collection of genomes**.
Given the genomes and a value $k$, it computes the $(\alpha^\ast, \beta_\mathrm{low}^\ast, \beta_\mathrm{high}^\ast)$ parameters for each genome's $k$-mer spectrum and outputs the results to a csv file.

## Building the tool

To build the tool, you need CMake (≥ 3.16), a C++ compiler (GCC ≥ 8 or Clang ≥ 5), and Python (≥ 3.10).

```bash
git clone --recursive https://github.com/medvedevgroup/CaPLa.git
cd CaPLa
mkdir build
cd build
cmake ..
make -j 8
```

## Usage

The tool [capla.sh](scripts/capla.sh) scans a given folder for genomes with `.fna` or `.fasta` extension and computes the CaPLa for each genome.
It takes as arguments the _path_ to the directory containing the genomes and an optional _k-mer size_ (default: 21).

For example, you can run the tool on the sample genomes already included in the `genomes` directory as follows:

```bash
./scripts/capla.sh genomes/ 21
```

The results will be saved in a file named `CaPLa.csv` inside the specified directory.
Each row of the file contains:
1. Genome name
2. Number of unique $k$-mers
3. $k$-mer length
4. $\alpha^\ast$
5. $\beta_\mathrm{low}^\ast$
6. $\beta_\mathrm{high}^\ast$

## Under the hood

The tool begins by preprocessing each genome: it uses the `process_fast` executable to remove non-`ACGT` characters and create a single string with a single header. It then builds a suffix array for each genome using the `mksary` executable from the `libdivsufsort` library.

Next, the tool computes the number of segments in the PLA of each genome's rank curve for each error bound $\varepsilon \in \{1, 2, 3, \dots, 1024\}$.
The rank curve captures the relationship between each $k$-mer and its rank in the suffix array.
This step is done using the `count_segments` executable, which implements O'Rourke's algorithm to find the PLA with the minimal number of segments.

Finally, the tool finds the CaPLa values for each genome using the `find_capla.py` script and writes the results to the `CaPLa.csv` file.

## Data used in paper submission
The data used in our paper submission is described in the [Reproducibility](Reproducibility) directory.
