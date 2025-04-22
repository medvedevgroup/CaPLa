BIN_PATH="$1"
MKSARY_PATH="$2"
SRC_PATH="$3"
KMER="$4"
OUTPUT_DIR="$(pwd)/"

count_processed=0
count_sa=0
count_seg=0
count_alpha_beta=0
# for file in *.fna; do
#     file_base="${file%.fna}"
#     file_base="${file_base%%.*}"
for file in *.fna *.fasta; do
    [[ -e "$file" ]] || continue  # skip if no matching files

    file_base=$(basename "$file")
    file_base="${file_base%%.*}"  # before first dot

    if [[ ! -f "$file_base.processed.fasta" ]]; then
        echo "Processing $file → ${file_base}.processed.fasta"
        "$BIN_PATH"process_fasta "$file" "$file_base"
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
        if [[ ! -f "$file.index.rourk.txt" ]]; then
            "$BIN_PATH"count_segments --genome_fasta="$file" --suffix_array="$file".sa -a -k "$KMER"
            count_seg=$((count_seg + 1))
        fi
    fi
done

FLAG=1
for file in *.index.rourk.txt; do
    file_base=$(basename "$file")
    file_base="${file_base%%.*}"  # before first dot
    # print echo for the first file and do not print for the rest
    if [[ $FLAG -eq 1 ]]; then
        python3 "${SRC_PATH}find_alpha_beta.py" "$file" "$file_base" X X X "$KMER" alpha_beta_list.csv F 
        FLAG=0
        echo $FLAG
    else
        python3 "${SRC_PATH}find_alpha_beta.py" "$file" "$file_base" X X X "$KMER" alpha_beta_list.csv T 
    fi
    count_alpha_beta=$((count_alpha_beta + 1))
done

echo "Processed files: $count_processed"
echo "Suffix arrays created: $count_sa"
echo "Segments-count files created: $count_seg"
echo "Number of genomes for which CaPLa is calculated: $count_alpha_beta"
