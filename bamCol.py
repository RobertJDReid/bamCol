#!/usr/bin/env python3
import argparse
import csv
import sys
from typing import Iterable, Tuple, List, Optional
import pysam

__version__ = "0.0.2"

def parse_pos_token(tok: str) -> Tuple[str, int]:
    """
    Parse CHR:POS (1-based) into (chrom, pos_int).
    """
    if ":" not in tok:
        raise ValueError(f"Position must be CHR:POS, got: {tok}")
    chrom, spos = tok.split(":", 1)
    pos = int(spos)
    if pos < 1:
        raise ValueError(f"Position must be 1-based positive integer, got: {pos}")
    return chrom, pos

def read_positions(pos_args: List[str], pos_file: Optional[str]) -> List[Tuple[str, int]]:
    positions = []
    for p in pos_args or []:
        positions.append(parse_pos_token(p))
    if pos_file:
        with open(pos_file) as fh:
            for line in fh:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                # tolerate TSV/CSV: chrom [tab|,] pos
                parts = [x.strip() for x in line.replace(",", "\t").split("\t") if x.strip()]
                if len(parts) < 2:
                    raise ValueError(f"Bad line in positions file: {line.strip()}")
                chrom = parts[0]
                pos = int(parts[1])
                if pos < 1:
                    raise ValueError(f"Position must be 1-based: {line.strip()}")
                positions.append((chrom, pos))
    if not positions:
        raise ValueError("No positions provided. Use --pos and/or --pos-file.")
    return positions

def pileup_one_position(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos1: int,
    min_mapq: int,
    min_bq: int,
    max_depth: int,
    ignore_overlaps: bool,
    ignore_orphans: bool,
    want_cigar: bool,
) -> Iterable[dict]:
    """
    Yield per-read records at a single 1-based reference position.
    """
    start0 = pos1 - 1
    try:
        cols = bam.pileup(
            chrom,
            start0,
            pos1,  # end is exclusive; pos1 ensures we cover exactly pos1
            truncate=True,
            stepper="samtools",
            min_mapping_quality=min_mapq,
            min_base_quality=0,  # we’ll filter per-base below so we can still report low-BQ if desired
            max_depth=max_depth,
            ignore_overlaps=ignore_overlaps,
            ignore_orphans=ignore_orphans,
        )
    except ValueError as e:
        # Often thrown if chromosome name not in BAM header
        raise ValueError(f"Chromosome {chrom!r} not found in BAM header? {e}")

    for col in cols:
        if col.reference_pos != start0:
            # With truncate=True this should already be true, but guard anyway.
            continue

        for pr in col.pileups:
            aln = pr.alignment

            # Skip secondary/supplementary/unmapped by default (pysam pileup already skips unmapped).
            if aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < min_mapq:
                continue

            is_del = pr.is_del
            is_refskip = pr.is_refskip  # e.g., spliced alignments in RNA-seq

            base = None
            base_qual = None
            qpos = pr.query_position  # None for deletions/refskips
            if not is_del and not is_refskip and qpos is not None:
                seq = aln.query_sequence
                quals = aln.query_qualities
                if seq is not None and 0 <= qpos < len(seq):
                    base = seq[qpos]
                if quals is not None and 0 <= qpos < len(quals):
                    base_qual = int(quals[qpos])

                # Filter by base quality if requested
                if base_qual is not None and base_qual < min_bq:
                    continue

            # Represent special cases
            call = base
            if is_del:
                call = "DEL"        # deletion w.r.t. reference base at this position
            elif is_refskip:
                call = "REFSKIP"    # e.g., 'N' CIGAR; no base at this ref pos

            rec = {
                "chrom": chrom,
                "pos": pos1,
                "read_id": aln.query_name,
                "call": call,
                "is_del": is_del,
                "is_refskip": is_refskip,
                "base_qual": base_qual if base_qual is not None else "",
                "mapq": aln.mapping_quality,
                "strand": "-" if aln.is_reverse else "+",
                "flag": aln.flag,
            }

            if want_cigar:
                rec["cigar"] = aln.cigarstring

            yield rec

def main():
    ap = argparse.ArgumentParser(
        description="Extract per-read allele calls at specific genome positions from a BAM using pysam."
    )
    ap.add_argument("bam", help="Input BAM (indexed .bai required)")
    ap.add_argument(
        "--pos",
        action="append",
        default=[],
        metavar="CHR:POS",
        help="1-based position like S288C_chrII:123927 (repeatable).",
    )
    ap.add_argument(
        "--pos-file",
        help="File of positions (TSV or CSV) with columns: chrom pos (1-based). Header and #comments allowed.",
    )
    ap.add_argument(
        "--out",
        default="-",
        help="Output CSV file (default: '-' for stdout).",
    )
    ap.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality (default: 0).")
    ap.add_argument("--min-bq", type=int, default=0, help="Minimum base quality (Phred) (default: 0).")
    ap.add_argument("--max-depth", type=int, default=100000, help="Max pileup depth (default: 100000).")
    ap.add_argument(
        "--no-ignore-overlaps",
        dest="ignore_overlaps",
        action="store_false",
        help="Count both mates if they overlap (default is to ignore overlap to avoid double counting).",
    )
    ap.add_argument(
        "--ignore-orphans",
        action="store_true",
        help="Ignore orphan reads (mates not properly paired).",
    )
    ap.add_argument(
        "--cigar",
        action="store_true",
        help="Include CIGAR string column in output.",
    )
    ap.add_argument(
        "--version", 
        action="version", 
        version=f"%(prog)s {__version__}"
    )

    args = ap.parse_args()

    positions = read_positions(args.pos, args.pos_file)

    # Open BAM
    try:
        bamf = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        ap.error(f"Failed to open BAM: {e}")

    # Verify index exists early
    if not bamf.has_index():
        ap.error("BAM index (.bai) not found or unreadable. Create with samtools index.")

    # Prepare output
    out_fh = sys.stdout if args.out == "-" else open(args.out, "w", newline="")

    # Build header depending on --cigar
    fieldnames = [
        "chrom",
        "pos",
        "read_id",
        "call",
        "is_del",
        "is_refskip",
        "base_qual",
        "mapq",
        "strand",
        "flag",
    ]
    if args.cigar:
        fieldnames.insert(fieldnames.index("flag"), "cigar")

    writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
    writer.writeheader()

    n_written = 0
    try:
        with bamf:
            for chrom, pos1 in positions:
                for rec in pileup_one_position(
                    bam=bamf,
                    chrom=chrom,
                    pos1=pos1,
                    min_mapq=args.min_mapq,
                    min_bq=args.min_bq,
                    max_depth=args.max_depth,
                    ignore_overlaps=args.ignore_overlaps,
                    ignore_orphans=args.ignore_orphans,
                    want_cigar=args.cigar,
                ):
                    writer.writerow(rec)
                    n_written += 1
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()

    # Optional: print a tiny summary to stderr
    print(f"Wrote {n_written} records", file=sys.stderr)

if __name__ == "__main__":
    main()
