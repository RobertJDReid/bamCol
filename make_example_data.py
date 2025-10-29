#!/usr/bin/env python3
import pysam
from pathlib import Path

example_dir = Path("example_data")
example_dir.mkdir(exist_ok=True)
bam_path = example_dir / "example.bam"
pos_path = example_dir / "positions.txt"

# Define a simple BAM header
header = {
    "HD": {"VN": "1.6"},
    "SQ": [{"LN": 1000, "SN": "S288C_chrI"}],
}

# Create a BAM file with a few example reads
with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read_001"
    a.query_sequence = "A" * 50
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 2930
    a.mapping_quality = 60
    a.cigar = ((0, 50),)
    a.query_qualities = pysam.qualitystring_to_array("I" * 50)
    outf.write(a)

    b = pysam.AlignedSegment()
    b.query_name = "read_002"
    b.query_sequence = "C" * 50
    b.flag = 147
    b.reference_id = 0
    b.reference_start = 2940
    b.mapping_quality = 59
    b.cigar = ((0, 50),)
    b.query_qualities = pysam.qualitystring_to_array("H" * 50)
    outf.write(b)

    c = pysam.AlignedSegment()
    c.query_name = "read_003"
    c.query_sequence = "G" * 50
    c.flag = 83
    c.reference_id = 0
    c.reference_start = 2945
    c.mapping_quality = 58
    c.cigar = ((0, 50),)
    c.query_qualities = pysam.qualitystring_to_array("J" * 50)
    outf.write(c)

# Index the BAM
pysam.index(str(bam_path))

# Create example positions file
pos_content = """# chrom    pos
S288C_chrI 2941
S288C_chrI 2947
"""
pos_path.write_text(pos_content)

print("Created example files in ./example_data/")
