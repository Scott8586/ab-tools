#!/usr/bin/env python3
"""
build_hcdr3_csv.py

Scans a directory of .fasta files (each containing ">H" / ">L" entries, as
produced by csv_to_fasta.py) and a directory of matching .pdb files (as
produced by run_abodybuilder2.sh), and writes a summary CSV describing the
location and sequence of each heavy chain's CDR3, using the `abnumber`
package for antibody numbering.

Usage:
    python build_hcdr3_csv.py <fasta_dir> <pdb_dir> <output_csv> [--scheme imgt]

Output CSV columns:
    index, fileid, pdb, ab_chains, chain, resi_start, resi_end, hcdr3
"""

import argparse
import csv
import os
import sys

from abnumber import Chain
from abnumber.exceptions import ChainParseError


def parse_fasta(path):
    """Parse a simple FASTA file into an OrderedDict of {header: sequence}."""
    records = {}
    header = None
    seq_chunks = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            records[header] = "".join(seq_chunks)

    return records


def get_heavy_chain_cdr3(vh_sequence, scheme):
    """Use abnumber to number the heavy chain sequence and pull out CDR3 info.

    Returns (resi_start, resi_end, cdr3_seq).
    """
    chain = Chain(vh_sequence, scheme=scheme)
    cdr3_dict = chain.cdr3_dict

    if not cdr3_dict:
        raise ValueError("abnumber found no CDR3 region for this sequence")

    positions = list(cdr3_dict.keys())
    resi_start = positions[0].number
    resi_end = positions[-1].number
    cdr3_seq = chain.cdr3_seq

    return resi_start, resi_end, cdr3_seq


def main():
    parser = argparse.ArgumentParser(
        description="Build a CSV summarizing heavy chain CDR3 location/sequence "
                    "from matching .fasta and .pdb directories."
    )
    parser.add_argument("fasta_dir", help="Directory containing <fileid>.fasta files")
    parser.add_argument("pdb_dir", help="Directory containing <fileid>.pdb files")
    parser.add_argument("output_csv", help="Path to write the output CSV")
    parser.add_argument(
        "--scheme",
        default="imgt",
        choices=["imgt", "chothia", "kabat", "aho"],
        help="Antibody numbering scheme to use with abnumber (default: imgt)",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.fasta_dir):
        print(f"Error: fasta directory '{args.fasta_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(args.pdb_dir):
        print(f"Error: pdb directory '{args.pdb_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    fasta_files = sorted(f for f in os.listdir(args.fasta_dir) if f.endswith(".fasta"))

    rows = []
    index = 0

    for fasta_filename in fasta_files:
        fileid = fasta_filename[: -len(".fasta")]
        fasta_path = os.path.join(args.fasta_dir, fasta_filename)
        pdb_path = os.path.join(args.pdb_dir, f"{fileid}.pdb")

        if not os.path.isfile(pdb_path):
            print(f"Skipping '{fileid}': no matching pdb file at '{pdb_path}'.")
            continue

        records = parse_fasta(fasta_path)
        vh_sequence = records.get("H")
        if not vh_sequence:
            print(f"Skipping '{fileid}': no '>H' entry found in '{fasta_path}'.")
            continue

        try:
            resi_start, resi_end, hcdr3 = get_heavy_chain_cdr3(vh_sequence, args.scheme)
        except (ChainParseError, ValueError) as e:
            print(f"Skipping '{fileid}': {e}")
            continue

        pdb_rel_path = os.path.relpath(pdb_path, start=args.pdb_dir)

        rows.append(
            {
                "index": index,
                "fileid": fileid,
                "pdb": pdb_rel_path,
                "ab_chains": "HL",
                "chain": "H",
                "resi_start": resi_start,
                "resi_end": resi_end,
                "hcdr3": hcdr3,
            }
        )
        index += 1

    fieldnames = ["index", "fileid", "pdb", "ab_chains", "chain", "resi_start", "resi_end", "hcdr3"]
    with open(args.output_csv, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Done. Wrote {len(rows)} row(s) to '{args.output_csv}'.")


if __name__ == "__main__":
    main()
