#!/usr/bin/env python3
"""
csv_to_fasta.py

Reads a CSV with columns: Id, Batch, Mutations, VH_VL
where VH_VL is "<VH_sequence>:<VL_sequence>".

For each row, writes a FASTA file named "<Id>.fasta" containing:
    >H
    <VH_sequence>
    >L
    <VL_sequence>

Usage:
    python csv_to_fasta.py input.csv [output_dir]

If output_dir is omitted, files are written to the current directory.
"""

import csv
import os
import sys


def sanitize_filename(name: str) -> str:
    """Strip characters that would be unsafe in a filename."""
    keep = "-_."
    return "".join(c for c in name if c.isalnum() or c in keep) or "unnamed"


def main():
    if len(sys.argv) < 2:
        print("Usage: python csv_to_fasta.py <input.csv> [output_dir]")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    os.makedirs(output_dir, exist_ok=True)

    written = 0
    skipped = 0

    with open(input_csv, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        required_cols = {"Id", "Batch", "Mutations", "VH_VL"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            print(f"Error: input CSV is missing required column(s): {', '.join(sorted(missing))}")
            sys.exit(1)

        for row_num, row in enumerate(reader, start=2):  # start=2: header is line 1
            rec_id = (row.get("Id") or "").strip()
            vh_vl = (row.get("VH_VL") or "").strip()

            if not rec_id:
                print(f"Row {row_num}: missing Id, skipping.")
                skipped += 1
                continue

            if ":" not in vh_vl:
                print(f"Row {row_num} (Id={rec_id}): VH_VL missing ':' separator, skipping.")
                skipped += 1
                continue

            vh_seq, vl_seq = vh_vl.split(":", 1)
            vh_seq = vh_seq.strip()
            vl_seq = vl_seq.strip()

            if not vh_seq or not vl_seq:
                print(f"Row {row_num} (Id={rec_id}): empty VH or VL sequence, skipping.")
                skipped += 1
                continue

            filename = f"{sanitize_filename(rec_id)}.fasta"
            filepath = os.path.join(output_dir, filename)

            with open(filepath, "w", encoding="utf-8") as out_f:
                out_f.write(">H\n")
                out_f.write(vh_seq + "\n")
                out_f.write(">L\n")
                out_f.write(vl_seq + "\n")

            written += 1

    print(f"Done. Wrote {written} FASTA file(s) to '{output_dir}'. Skipped {skipped} row(s).")


if __name__ == "__main__":
    main()
