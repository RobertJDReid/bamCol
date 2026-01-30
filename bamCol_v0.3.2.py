#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
import warnings
from typing import Iterable, Tuple, List, Optional
import pysam
from multiprocessing import Pool

__version__ = "0.3.2"


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


def parse_region_token(token: str) -> Tuple[str, int, int]:
    """
    Parse 'CHR:START-END' or 'CHR:START..END' into (chrom, start, end), inclusive, 1-based.
    """
    # Support both - and .. as separators
    t = token.strip().replace("..", "-")
    
    # Pattern: CHR:START-END
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", t)
    if not m:
        raise ValueError(f"Bad region '{token}'. Expected format: CHR:START-END (e.g., chr1:1000-2000)")
    
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    
    # Auto-swap if start > end
    if start > end:
        start, end = end, start
    
    # Validate 1-based positions
    if start < 1 or end < 1:
        raise ValueError(f"Region positions must be >= 1: {token}")
    
    return chrom, start, end


def in_regions(chrom: str, pos: int, regions: List[Tuple[str, int, int]]) -> bool:
    """
    Check if (chrom, pos) falls within ANY region in the regions list.
    Returns True if position is in any region, False otherwise.
    """
    for rchrom, start, end in regions:
        if chrom == rchrom and start <= pos <= end:
            return True
    return False


def read_vcf_positions(vcf_path: str) -> List[Tuple[str, int]]:
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
        positions: List[Tuple[str, int]] = []

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

            # Validate position (defensive programming)
            if rec.pos < 1:
                warnings.warn(f"Skipping invalid position in VCF: {rec.contig}:{rec.pos}")
                continue

            positions.append((rec.contig, rec.pos))

        vcf.close()
        return positions
        
    except Exception as pysam_error:
        # Fall back to manual parsing for regular gzip files
        # (pysam fails on regular gzip, only works with bgzip)
        if not vcf_path.endswith('.gz'):
            # If it's not gzipped and pysam failed, re-raise the original error
            raise pysam_error
            
        positions: List[Tuple[str, int]] = []
        with gzip.open(vcf_path, 'rt') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    warnings.warn(f"Skipping malformed VCF line (< 7 fields): {line.strip()[:100]}")
                    continue
                    
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    warnings.warn(f"Skipping VCF line with invalid position: {line.strip()[:100]}")
                    continue
                
                ref = parts[3]
                alt = parts[4]
                filt = parts[6]
                
                # Validate position
                if pos < 1:
                    warnings.warn(f"Skipping invalid position in VCF: {chrom}:{pos}")
                    continue
                
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


def read_positions(
    pos_args: List[str],
    pos_file: Optional[str],
    vcf_file: Optional[str],
    include_regions: List[Tuple[str, int, int]],
    exclude_regions: List[Tuple[str, int, int]]
) -> List[Tuple[str, int]]:
    """
    Collect positions from --pos arguments, --pos-file, and --vcf-file.
    Apply include/exclude region filtering.
    Returns deduplicated list of (chrom, pos) tuples.
    """
    positions = []
    
    # From --pos arguments
    for p in pos_args or []:
        positions.append(parse_pos_token(p))
    
    # From --pos-file
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
                try:
                    pos = int(parts[1])
                except ValueError as e:
                    raise ValueError(f"Invalid position value in line: {line.strip()}") from e
                if pos < 1:
                    raise ValueError(f"Position must be 1-based: {line.strip()}")
                positions.append((chrom, pos))
    
    # From --vcf-file
    if vcf_file:
        vcf_positions = read_vcf_positions(vcf_file)
        if not vcf_positions:
            warnings.warn(
                f"No valid positions found in VCF file '{vcf_file}'. "
                f"Only biallelic SNPs with FILTER=PASS are used."
            )
        positions.extend(vcf_positions)
    
    if not positions:
        raise ValueError("No positions provided. Use --pos, --pos-file, and/or --vcf-file.")
    
    # Apply region filtering
    filtered = []
    for chrom, pos in positions:
        # If include_regions specified, position must be in at least one
        if include_regions and not in_regions(chrom, pos, include_regions):
            continue
        # If exclude_regions specified, position must not be in any
        if exclude_regions and in_regions(chrom, pos, exclude_regions):
            continue
        filtered.append((chrom, pos))
    
    # Deduplicate and sort
    filtered = sorted(set(filtered), key=lambda x: (x[0], x[1]))
    
    return filtered


def pileup_one_position(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos1: int,
    min_mapq: int,
    min_bq: int,
    max_depth: int,
    ignore_overlaps: bool,
    ignore_orphans: bool,
    include_secondary: bool,
    include_supplementary: bool,
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
            min_base_quality=0,  # we'll filter per-base below so we can still report low-BQ if desired
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

            # Skip secondary/supplementary alignments unless requested
            if aln.is_secondary and not include_secondary:
                continue
            if aln.is_supplementary and not include_supplementary:
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
                "read_pos": qpos if qpos is not None else "",
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


def worker_pileup(args) -> List[dict]:
    """
    Worker for a single position.

    Opens its own AlignmentFile (needed for safe multiprocessing),
    runs pileup_one_position, and returns a list of records.
    """
    (
        bam_path,
        chrom,
        pos1,
        min_mapq,
        min_bq,
        max_depth,
        ignore_overlaps,
        ignore_orphans,
        include_secondary,
        include_supplementary,
        want_cigar,
    ) = args

    records: List[dict] = []
    bamf = pysam.AlignmentFile(bam_path, "rb")
    try:
        for rec in pileup_one_position(
            bam=bamf,
            chrom=chrom,
            pos1=pos1,
            min_mapq=min_mapq,
            min_bq=min_bq,
            max_depth=max_depth,
            ignore_overlaps=ignore_overlaps,
            ignore_orphans=ignore_orphans,
            include_secondary=include_secondary,
            include_supplementary=include_supplementary,
            want_cigar=want_cigar,
        ):
            records.append(rec)
    finally:
        bamf.close()
    return records


def main():
    ap = argparse.ArgumentParser(
        description="Extract per-read allele calls at specific genome positions from a BAM using pysam.",
        epilog="Note: For complex region filtering (>10 regions), consider pre-filtering "
               "your VCF with bcftools for better performance."
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
        "--vcf-file",
        help="Input VCF file. Biallelic SNP positions with FILTER=PASS will be used as query sites (VCF POS is 1-based).",
    )
    ap.add_argument(
        "--include-region",
        action="append",
        default=[],
        metavar="CHR:START-END",
        help="Include only positions within region (repeatable). Format: CHR:START-END, e.g., chr1:1000-2000. "
             "Recommended for up to ~10 regions.",
    )
    ap.add_argument(
        "--exclude-region",
        action="append",
        default=[],
        metavar="CHR:START-END",
        help="Exclude positions within region (repeatable). Format: CHR:START-END, e.g., chr1:1000-2000. "
             "Recommended for up to ~10 regions.",
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
        "--include-secondary",
        action="store_true",
        help="Include secondary alignments (FLAG 0x100) in pileup output.",
    )
    ap.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Include supplementary alignments (FLAG 0x800) in pileup output.",
    )
    ap.add_argument(
        "--cigar",
        action="store_true",
        help="Include CIGAR string column in output.",
    )
    ap.add_argument(
        "--process",
        type=int,
        default=None,
        help="Number of worker processes (default: use all available CPUs).",
    )
    ap.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = ap.parse_args()

    # Parse include/exclude regions
    include_regions = []
    exclude_regions = []
    
    try:
        for region_str in args.include_region:
            include_regions.append(parse_region_token(region_str))
        
        for region_str in args.exclude_region:
            exclude_regions.append(parse_region_token(region_str))
    except ValueError as e:
        ap.error(str(e))

    # Warn if using many regions (performance consideration)
    total_regions = len(include_regions) + len(exclude_regions)
    if total_regions > 10:
        warnings.warn(
            f"Using {total_regions} regions may impact performance. "
            f"For large region lists, consider pre-filtering with bcftools."
        )

    # Resolve positions from all sources with region filtering
    try:
        positions = read_positions(
            args.pos,
            args.pos_file,
            args.vcf_file,
            include_regions,
            exclude_regions
        )
    except Exception as e:
        ap.error(str(e))

    # Verify BAM + index once in the main process
    try:
        bam_test = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        ap.error(f"Failed to open BAM: {e}")

    if not bam_test.has_index():
        bam_test.close()
        ap.error("BAM index (.bai) not found or unreadable. Create with samtools index.")
    bam_test.close()

    # Determine number of processes
    cpu_count = os.cpu_count() or 1
    if args.process is None:
        n_proc = cpu_count
    else:
        n_proc = args.process
    if n_proc < 1:
        ap.error(f"--process must be >= 1 (got {n_proc})")

    print(
        f"Using {n_proc} worker process(es) out of {cpu_count} available CPU(s).",
        file=sys.stderr,
    )
    if include_regions or exclude_regions:
        print(
            f"Region filtering: {len(include_regions)} include, {len(exclude_regions)} exclude",
            file=sys.stderr,
        )
    print(
        f"Processing {len(positions)} position(s).",
        file=sys.stderr,
    )

    # Prepare output
    try:
        out_fh = sys.stdout if args.out == "-" else open(args.out, "w", newline="")

        # Build header depending on --cigar
        fieldnames = [
            "chrom",
            "pos",
            "read_id",
            "read_pos",
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

        # Single-process path: behave like the original code
        if n_proc == 1:
            bamf = pysam.AlignmentFile(args.bam, "rb")
            try:
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
                        include_secondary=args.include_secondary,
                        include_supplementary=args.include_supplementary,
                        want_cigar=args.cigar,
                    ):
                        writer.writerow(rec)
                        n_written += 1
            finally:
                bamf.close()

        # Multi-process path
        else:
            # Build task list: one job per position
            task_args = [
                (
                    args.bam,
                    chrom,
                    pos1,
                    args.min_mapq,
                    args.min_bq,
                    args.max_depth,
                    args.ignore_overlaps,
                    args.ignore_orphans,
                    args.include_secondary,
                    args.include_supplementary,
                    args.cigar,
                )
                for chrom, pos1 in positions
            ]

            # Use imap to preserve order of positions in output
            with Pool(processes=n_proc) as pool:
                for records in pool.imap(worker_pileup, task_args):
                    for rec in records:
                        writer.writerow(rec)
                        n_written += 1

        # Tiny summary to stderr
        print(f"Wrote {n_written} records", file=sys.stderr)

    finally:
        if out_fh is not sys.stdout:
            out_fh.close()


if __name__ == "__main__":
    main()
