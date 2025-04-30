SCRIPT_PATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
BIN_PATH="$SCRIPT_PATH/../build"
MKSARY_PATH="$BIN_PATH/mksary"
SRC_PATH="$SCRIPT_PATH/../src"
KMER="${2:-21}"

if [[ -z "$1" ]]; then
    echo "Usage: $0 <directory> [kmer]"
    echo "  directory: Directory containing the .fna or .fasta files."
    echo "  kmer: Optional k-mer size (default is 21)."
    exit 1
fi

if [[ -d "$1" ]]; then
    cd "$1" || exit 1
else
    echo "Directory $1 does not exist."
    exit 1
fi

if [[ ! -d "$BIN_PATH/venv" ]]; then
    echo "Creating virtual environment in $BIN_PATH/venv"
    python3 -m venv "$BIN_PATH/venv"
    source "$BIN_PATH/venv/bin/activate"
    python3 -m pip install "pandas>=2.2.3" "scipy>=1.15.2"
else
    source "$BIN_PATH/venv/bin/activate"
fi

count_processed=0
count_sa=0
count_seg=0
count_capla=0

shopt -s nullglob # Do not match empty glob patterns

for file in *.fna *.fasta; do
    [[ -e "$file" ]] || continue  # skip if no matching files

    file_base=$(basename "$file")
    file_base="${file_base%%.*}"  # before first dot

    if [[ ! -f "$file_base.processed.fasta" ]]; then
        echo "Processing $file → ${file_base}.processed.fasta"
        "$BIN_PATH"/process_fasta "$file" "$file_base"
        count_processed=$((count_processed + 1))
    fi
done

for file in *.processed.fasta; do
    bases_count=$(tail -n +2 "$file" | wc -c)
    if (( bases_count >= 10000 )); then
        if [[ ! -f "$file.sa" ]]; then
            "$MKSARY_PATH" "$file" "$file.sa"
            count_sa=$((count_sa + 1))
        fi
        if [[ ! -f "$file.segments.txt" ]]; then
            "$BIN_PATH"/count_segments --genome_fasta="$file" --suffix_array="$file".sa -a -k "$KMER"
            count_seg=$((count_seg + 1))
        fi
    fi
done

FLAG=1
for file in *.segments.txt; do
    file_base=$(basename "$file")
    file_base="${file_base%%.*}"  # before first dot
    # print echo for the first file and do not print for the rest
    if [[ $FLAG -eq 1 ]]; then
        python3 "${SRC_PATH}/find_capla.py" "$file" "$file_base" "$KMER" CaPLa.csv F 
        FLAG=0
        echo $FLAG
    else
        python3 "${SRC_PATH}/find_capla.py" "$file" "$file_base" "$KMER" CaPLa.csv T 
    fi
    count_capla=$((count_capla + 1))
done

echo "Processed files: $count_processed"
echo "Suffix arrays created: $count_sa"
echo "Segments-count files created: $count_seg"
echo "Number of genomes for which CaPLa is calculated: $count_capla"
