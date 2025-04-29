# CaPLa

Canonical piecewise linear approximation (CaPLa) represents the tightest bound of two power-laws on average segment length for applications that use piecewise linear approximation.
Here, we provide a script that can be used to find CaPLa for a set of genomes.
Essentially, for each genome, we find $(\alpha^*, \beta_{low}^*, \beta_{high}^*)$ and write them in a csv file.

## Installation

```bash
git clone --recursive https://github.com/medvedevgroup/CaPLa.git
cd CaPLa
mkdir build
cd build
cmake ..
make -j 8
```

## Usage

You can use the script `find_CaPLa.sh` to find CaPLa for a set of genomes.
We assume all the genome files have either `.fna` or `.fasta` extension.
The script will look for all the files in the current directory with these extensions and will use them as input genomes.
The script takes four arguments:

1. The path to the build directory.
2. The path to the `mksary` executable.
3. The path to the source directory.
4. $k$-mer length

At first, the script will process each of the genomes to remove any non `ACGT` characters and concatenate them into a single string with a single header.
Then, it will create suffix arrays for each genome using the `mksary` executable.
The `mksary` executable is a part of the `libdivsufsort` library, which is used to create suffix arrays.
Afterwards, it will count the number of segments for each genome using the `count_segments` executable for $\epsilon \in \{1, 2, 3, \dots, 1024\}$.
The `count_segments` executable uses O'Rourke's algorithm to find the optimal number of segments of a genome's rank curve at the corresponding $\epsilon$.
The rank curve captures the relationship between each $k$-mer and it's rank in the suffix array.
Finally, it will find the CaPLa values for each genome using the `find_CaPLa.py` file, and write the results to a csv file named `CaPLa.csv`.

Let's say you have a set of genomes in the folder `genomes/` and you want to find CaPLa for all of them. You can do this by running the following command:

```bash
cd genomes
../scripts/capla.sh ../build ../build/mksary ../src 21 
```

The results will be written to a file called `CaPLa.csv` inside the `genomes` directory.
This file will contain the following columns:

1. Genome name
2. Number of unique $k$-mers
3. $k$-mer length
4. $\alpha^*$
5. $\beta_{low}^*$
6. $\beta_{high}^*$
