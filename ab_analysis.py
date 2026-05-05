#!/usr/bin/env python3
"""
ab_analysis.py — Analyse antibody sequences from a FASTA file.

Author:     Scott Presnell <srp@presnellgroup.net>
Date:       2026-05-04
Version:    0.9.0
License:    MIT

Usage:
    python ab_analysis.py sequences.fasta
    python ab_analysis.py --scheme kabat sequences.fasta 
    python ab_analysis.py --scheme imgt --csv sequences.fasta 

Requirements:
    conda install -c bioconda abnumber
"""

import argparse
import sys

try:
    from abnumber import Chain
    from abnumber.exceptions import ChainParseError
except ImportError:
    sys.exit(
        "ERROR: abnumber is not installed.\n"
        "Install it with:  conda install -c bioconda abnumber"
    )


# ── FASTA parser (no Biopython dependency) ────────────────────────────────────

def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:], []
            elif line:
                parts.append(line.upper())
    if header is not None:
        yield header, "".join(parts)


# ── CDR extraction ─────────────────────────────────────────────────────────────

REGION_NAMES = ["fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3"]

def analyse_sequence(name, seq, scheme):
    """
    Return a dict with chain type and CDR info, or an error string.
    """
    try:
        chain = Chain(seq, scheme=scheme, assign_germline=True)
    except ChainParseError as e:
        return {"error": str(e)}
    except Exception as e:
        return {"error": f"Unexpected error: {e}"}

    regions = {}
    for region in REGION_NAMES:
        region_seq = getattr(chain, f"{region}_seq", "") or ""
        regions[region] = {"seq": region_seq, "length": len(region_seq)}

    return {
        "length_overall": str(len(seq)),
        "species": getattr(chain, "species"),
        "germline": getattr(chain, "v_gene"),
        "chain_type": chain.chain_type,   # 'H', 'K', or 'L'
        "regions": regions,
    }


# ── Output formatters ──────────────────────────────────────────────────────────

def print_human_readable(results, scheme):
    print(f"\n{'='*60}")
    print(f"  CDR Analysis   |   Numbering scheme: {scheme.upper()}")
    print(f"{'='*60}\n")

    for name, data in results:
        print(f"► {name}")
        if "error" in data:
            print(f"    Could not parse: {data['error']}\n")
            continue

        chain_label = {"H": "Heavy", "K": "Kappa (light)", "L": "Lambda (light)"}.get(
            data["chain_type"], data["chain_type"]
        )
        germline = data["germline"]
        species = data["species"]
        print(f"  Chain type : {chain_label}")
        print(f"  Species    : {species}")
        print(f"  Germline   : {germline}")
        print(f"  {'CDR':<8} {'Sequence':<40} {'Length':>6}")
        print(f"  {'-'*8} {'-'*40} {'-'*6}")
        for region, info in data["regions"].items():
            seq_display = info["seq"] if info["seq"] else "(not found)"
            print(f"  {region.upper():<8} {seq_display:<40} {info['length']:>6}")
        print()


def print_csv(results, scheme):
    header = ",".join(["sequence_id", "chain_type", "species", "germline", "scheme",
                        "fr1_seq", "fr1_length",
                        "cdr1_seq", "cdr1_length",
                        "fr2_seq", "fr2_length",
                        "cdr2_seq", "cdr2_length",
                        "fr3_seq", "fr3_length",
                        "cdr3_seq", "cdr3_length",
                        "length_overall",
                        "error"])
    print(header)

    for name, data in results:
        if "error" in data:
            row = [name, "", scheme, "", "", "", "", "", "", "", "", "", "", "", "", data["error"]]
        else:
            row = [name, data["chain_type"], data["species"], data["germline"], scheme]
            for region in REGION_NAMES:
                info = data["regions"][region]
                row += [info["seq"], str(info["length"])]
            row.append(data["length_overall"])
            row.append("")
        print(",".join(row))


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Analyse antibody sequences from a FASTA file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Author: Scott Presnell <srp@presnellgroup.net>"
    )
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument(
        "--scheme",
        default="imgt",
        choices=["imgt", "kabat", "chothia", "martin", "aho"],
        help="Numbering scheme (default: imgt)",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Output as comma separated values instead of human-readable text",
    )
    args = parser.parse_args()

    # Parse and analyse
    try:
        records = list(parse_fasta(args.input_file))
    except FileNotFoundError:
        sys.exit(f"ERROR: File not found: {args.input_file}")

    if not records:
        sys.exit(f"ERROR: No sequences found in {args.input_file}")

    results = []
    for name, seq in records:
        results.append((name, analyse_sequence(name, seq, args.scheme)))

    # Print summary to stderr so it doesn't pollute TSV output
    n_ok = sum(1 for _, d in results if "error" not in d)
    print(
        f"Processed {len(results)} sequence(s): {n_ok} OK, "
        f"{len(results) - n_ok} failed",
        file=sys.stderr,
    )

    if args.csv:
        print_csv(results, args.scheme)
    else:
        print_human_readable(results, args.scheme)


if __name__ == "__main__":
    main()
