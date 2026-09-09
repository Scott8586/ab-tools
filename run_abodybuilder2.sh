#!/usr/bin/env bash
#
# run_abodybuilder2.sh
#
# Runs ABodyBuilder2 on every .fasta file in an input directory, writing
# predicted structures to an output directory. Skips any file whose
# corresponding .pdb already exists in the output directory.
#
# Usage:
#   ./run_abodybuilder2.sh <input_dir> <output_dir>

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>" >&2
    exit 1
fi

input_dir="$1"
output_dir="$2"

if [ ! -d "$input_dir" ]; then
    echo "Error: input directory '$input_dir' does not exist." >&2
    exit 1
fi

mkdir -p "$output_dir"

shopt -s nullglob
fasta_files=("$input_dir"/*.fasta)
shopt -u nullglob

if [ "${#fasta_files[@]}" -eq 0 ]; then
    echo "No .fasta files found in '$input_dir'."
    exit 0
fi

processed=0
skipped=0
failed=0

for fasta_path in "${fasta_files[@]}"; do
    filename="$(basename "$fasta_path")"
    fileid="${filename%.fasta}"
    pdb_path="$output_dir/${fileid}.pdb"

    if [ -f "$pdb_path" ]; then
        echo "Skipping '$fileid': '$pdb_path' already exists."
        skipped=$((skipped + 1))
        continue
    fi

    echo "Running ABodyBuilder2 on '$fileid'..."
    if Abodybuilder2 -f "$input_dir/${fileid}.fasta" -o "$pdb_path"; then
        processed=$((processed + 1))
    else
        echo "Error: ABodyBuilder2 failed for '$fileid'." >&2
        failed=$((failed + 1))
    fi
done

echo "Done. Processed: $processed, Skipped: $skipped, Failed: $failed."

if [ "$failed" -gt 0 ]; then
    exit 1
fi
