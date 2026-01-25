#!/usr/bin/env python3
import argparse
import csv
import sys
from typing import Iterable, Tuple, List, Optional
import pysam
##
import re

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
            min_base_quality=0,  # weâ€™ll filter per-base below so we can still report low-BQ if desired
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
                call = "DEL"   # deletion w.r.t. reference base at this position
            elif is_refskip:
                call = "REFSKIP"  # e.g., 'N' cigar; no base at this ref pos

            yield {
                "chrom": chrom,
                "pos": pos1,
                "read_id": aln.query_name,
                "call": call,
                "is_del": is_del,
                "is_refskip": is_refskip,
                "base_qual": base_qual if base_qual is not None else "",
                "mapq": aln.mapping_quality,
                "strand": "-" if aln.is_reverse else "+",
                "cigar": aln.cigarstring,
                "flag": aln.flag,
            }


##
_REGION_RE = re.compile(r"^([^:]+):(\d+)(?:-|\.\.)(\d+)$")


def parse_region_token(token: str) -> tuple[str, int, int]:
    """
    Parse 'CHR:START-END' or 'CHR:START..END' into (chrom, start, end), inclusive, 1-based.
    """
    t = token.strip().replace("..", "-")
    m = _REGION_RE.match(t)
    if not m:
        raise ValueError(f"Bad region '{token}'. Expected CHR:START-END")
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    if start > end:
        start, end = end, start
    return chrom, start, end


def read_region_file(path: str) -> list[tuple[str, int, int]]:
    """
    Reads regions from a TSV/CSV/whitespace-like file with columns: chrom start end (1-based, inclusive).
    Header allowed; lines beginning with # are ignored.
    Delimiter is auto-handled in a simple way:
      - if line contains '\t' -> split by tab
      - elif contains ',' -> split by comma
      - else -> split by whitespace
    """
    regions: list[tuple[str, int, int]] = []
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue

            # skip header-ish lines
            low = s.lower()
            if "chrom" in low and ("start" in low or "pos" in low) and "end" in low:
                continue

            if "\t" in s:
                parts = s.split("\t")
            elif "," in s:
                parts = s.split(",")
            else:
                parts = re.split(r"\s+", s)

            if len(parts) < 3:
                continue

            chrom = parts[0].strip()
            start = int(parts[1])
            end = int(parts[2])
            if start > end:
                start, end = end, start
            regions.append((chrom, start, end))

    return regions


def in_regions(chrom: str, pos: int, regions: list[tuple[str, int, int]]) -> bool:
    """
    True if (chrom,pos) falls within ANY region in regions.
    """
    for rchrom, start, end in regions:
        if chrom == rchrom and start <= pos <= end:
            return True
    return False


def read_vcf_positions(vcf_path: str) -> list[tuple[str, int]]:
    """
    Read positions from VCF file (handles .vcf, .vcf.gz with bgzip or regular gzip).
    Minimal default protocol:
      - use only biallelic SNPs
      - require FILTER=PASS (or empty filter)  [common VCF convention]
    Returns (contig, pos) where pos is 1-based (VCF standard).
    """
    import gzip
    
    # Try pysam first (works with bgzip and uncompressed VCF)
    try:
        vcf = pysam.VariantFile(vcf_path)
        positions: list[tuple[str, int]] = []

        # If not indexed, fetch() may fail; iterate directly.
        try:
            it = vcf.fetch()
        except Exception:
            it = vcf

        for rec in it:
            # PASS-only (handle writers that represent PASS as empty filter)
            fkeys = set(rec.filter.keys())
            if fkeys and "PASS" not in fkeys:
                continue

            # biallelic SNPs only
            if rec.alts is None or len(rec.alts) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue

            positions.append((rec.contig, int(rec.pos)))

        vcf.close()
        return positions
        
    except Exception as pysam_error:
        # Fall back to manual parsing for regular gzip files
        # (pysam fails on regular gzip, only works with bgzip)
        if not vcf_path.endswith('.gz'):
            # If it's not gzipped and pysam failed, re-raise the original error
            raise pysam_error
            
        positions: list[tuple[str, int]] = []
        with gzip.open(vcf_path, 'rt') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue
                    
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                filt = parts[6]
                
                # PASS filter only (VCF uses '.' for missing filter, treat as PASS)
                if filt != 'PASS' and filt != '.':
                    continue
                    
                # Biallelic SNPs only (skip if multiple alts separated by comma)
                if ',' in alt:
                    continue
                if len(ref) != 1 or len(alt) != 1:
                    continue
                    
                positions.append((chrom, pos))
        
        return positions


def collect_positions_with_regions(args) -> list[tuple[str, int]]:
    """
    Combine --pos, --pos-file, and --vcf-file into a single position list.
    Then apply include/exclude region gates.
    """
    positions: list[tuple[str, int]] = []

    # existing --pos
    for token in args.pos:
        chrom, pos = token.split(":", 1)
        positions.append((chrom, int(pos)))

    # existing --pos-file (keep using your current implementation if you already have one)
    if args.pos_file:
        # If you already have a function for this, call it instead.
        # This is a minimal inline reader similar to region-file handling.
        with open(args.pos_file, "r", encoding="utf-8") as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith("#"):
                    continue
                low = s.lower()
                if "chrom" in low and ("pos" in low or "position" in low):
                    continue

                if "\t" in s:
                    parts = s.split("\t")
                elif "," in s:
                    parts = s.split(",")
                else:
                    parts = re.split(r"\s+", s)

                if len(parts) < 2:
                    continue
                chrom = parts[0].strip()
                pos = int(parts[1])
                positions.append((chrom, pos))

    # NEW --vcf-file
    if args.vcf_file:
        positions.extend(read_vcf_positions(args.vcf_file))

    if not positions:
        raise SystemExit("No positions provided. Use --pos / --pos-file / --vcf-file.")

    # NEW include/exclude regions
    include_regions: list[tuple[str, int, int]] = []
    exclude_regions: list[tuple[str, int, int]] = []

    # --include-region (repeatable)
    for token in args.include_region:
        include_regions.append(parse_region_token(token))

    # --include-region-file
    if args.include_region_file:
        include_regions.extend(read_region_file(args.include_region_file))

    # --exclude-region (repeatable)
    for token in args.exclude_region:
        exclude_regions.append(parse_region_token(token))

    # Apply gating:
    # 1) if include_regions provided, keep only positions inside them
    # 2) always drop positions inside exclude_regions
    filtered: list[tuple[str, int]] = []
    for chrom, pos in positions:
        if include_regions and not in_regions(chrom, pos, include_regions):
            continue
        if exclude_regions and in_regions(chrom, pos, exclude_regions):
            continue
        filtered.append((chrom, pos))

    # de-duplicate + stable order
    filtered = sorted(set(filtered), key=lambda x: (x[0], x[1]))
    return filtered

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
    ##
    ap.add_argument(
        "--vcf-file",
        default=None,
        help="Input VCF file. SNP positions will be used as query sites (VCF POS is 1-based).",
    )
    ap.add_argument(
        "--include-region",
        action="append",
        default=[],
        metavar="CHR:START-END",
        help="Include only positions within region (repeatable). Format: CHR:START-END (1-based, inclusive).",
    )
    ap.add_argument(
        "--include-region-file",
        default=None,
        help="File of regions with columns: chrom start end (1-based, inclusive). Header and #comments allowed.",
    )
    ap.add_argument(
        "--exclude-region",
        action="append",
        default=[],
        metavar="CHR:START-END",
        help="Exclude positions within region (repeatable). Format: CHR:START-END (1-based, inclusive).",
    )


    args = ap.parse_args()

    ##
    #positions = read_positions(args.pos, args.pos_file)
    positions = collect_positions_with_regions(args)


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
    writer = csv.DictWriter(
        out_fh,
        fieldnames=[
            "chrom",
            "pos",
            "read_id",
            "call",
            "is_del",
            "is_refskip",
            "base_qual",
            "mapq",
            "strand",
            "cigar",
            "flag",
        ],
    )
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
