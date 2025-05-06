SCRIPT_PATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
BIN_PATH="$SCRIPT_PATH/../build"
MKSARY_PATH="$BIN_PATH/mksary"
SRC_PATH="$SCRIPT_PATH/../src"

# KMER="${2:-21}"
# EPS_MAX="${3:-1024}"
# USE_ONLY_TWO_POWER_EPS="${4:-no}"

KMER="21"
EPS_MAX="1024"
USE_ONLY_TWO_POWER_EPS="no"
directory=""

usage() {                                      # Function: Print a help message.
    echo "CaPLa Usage: -d directory [ -k kmer_size ] [ -e eps_max ] [ -t ] " 

    # echo "Usage: $0 <directory> [kmer]"
    echo "  -d directory: Directory containing the .fna or .fasta files."
    echo "  -k kmer_size: [INT] Optional k-mer size (default is 21)."
    echo "  -e eps_max: [INT] Optional maximum epsilon value (default is 1024)."
    echo "  -t: Optional flag to make CaPLa use only the epsilon values that are power of 2 instead of consecutive ones. If not specified, it will use consecutive values upto maximum epsilon."
    echo "  -h: Print help and exit."
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}

while getopts ":d:k:e:th" opt; do
  case $opt in
    d) directory="$OPTARG" ;;
    k) KMER="$OPTARG" ;;
    e) EPS_MAX="$OPTARG" ;;
    t) USE_ONLY_TWO_POWER_EPS="yes" ;;
    h) usage; exit 0 ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit_abnormal ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit_abnormal ;;
  esac
done

# shift $((OPTIND -1)) # Shift off the options and optional --.

if [[ -z "$directory" ]]; then
    echo "Error: Directory not specified."
    exit_abnormal
fi
if [[ ! -d "$directory" ]]; then
    echo "Error: Directory $directory does not exist."
    exit_abnormal
fi

set -e # Exit immediately if a command exits with a non-zero status

# if [[ -z "$1" ]]; then
#     exit 1
# fi

cd "$directory"

if [[ ! -d "$BIN_PATH/venv" ]]; then
    echo "Creating virtual environment in $BIN_PATH/venv"
    python3 -m venv "$BIN_PATH/venv"
    source "$BIN_PATH/venv/bin/activate"
    python3 -m pip install "pandas>=2.2.3" "scipy>=1.15.2"
else
    source "$BIN_PATH/venv/bin/activate"
fi

# EPS_FLAG='-a'

# if [[ $USE_ONLY_TWO_POWER_EPS == "yes" ]]; then
#     EPS_FLAG=''
# fi


count_input_files=0
count_processed=0
count_sa=0
count_seg=0
count_capla=0

shopt -s nullglob # Do not match empty glob patterns

for file in *.fna *.fasta; do
    [[ -e "$file" ]] || continue  # skip if no matching files

    file_base=$(basename "$file")
    file_base="${file_base%%.*}"  # before first dot
    if [[ "$file" != "$file_base.processed.fasta" ]]; then
        count_input_files=$((count_input_files + 1))        
    fi
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
        if [[ ! -f "$file."$KMER".segments.txt" ]]; then
            # echo "$BIN_PATH"/count_segments --genome_fasta="$file" --suffix_array="$file".sa "$EPS_FLAG" -k "$KMER" -e "$EPS_MAX" 
            if [[ $USE_ONLY_TWO_POWER_EPS == "yes" ]]; then
                "$BIN_PATH"/count_segments --genome_fasta="$file" --suffix_array="$file".sa  -k "$KMER" -e "$EPS_MAX"
            else
                "$BIN_PATH"/count_segments --genome_fasta="$file" --suffix_array="$file".sa -a -k "$KMER" -e "$EPS_MAX"
            fi
            count_seg=$((count_seg + 1))
        fi
    fi
done

for file in *$KMER.segments.txt; do
    file_base=$(basename "$file")
    file_base="${file_base%%.*}"  # before first dot
    # print echo for the first file and do not print for the rest
    python3 "${SRC_PATH}/find_capla.py" "$file" "$file_base" "$KMER" CaPLa.csv 
    count_capla=$((count_capla + 1))
done

echo "Input files: $count_input_files"
echo "Processed files: $count_processed"
echo "Suffix arrays created: $count_sa"
echo "Segments-count files created: $count_seg"
echo "Number of genomes for which CaPLa is calculated: $count_capla"
