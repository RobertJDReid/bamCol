#!/usr/bin/env python3
import pysam
from pathlib import Path

example_dir = Path("example_data")
example_dir.mkdir(exist_ok=True)

# File paths
bam_path = example_dir / "example.bam"
pos_path = example_dir / "positions.txt"
vcf_path = example_dir / "example.vcf"
malformed_vcf_path = example_dir / "malformed.vcf"

print("Creating example data files...")

# ============================================================================
# 1. Create BAM file with example reads
# ============================================================================
print("  Creating BAM file...")

# Define a simple BAM header
header = {
    "HD": {"VN": "1.6"},
    "SQ": [
        {"LN": 10000, "SN": "S288C_chrI"},
        {"LN": 10000, "SN": "S288C_chrII"}
    ],
}

# Create a BAM file with reads covering various positions
with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
    # Read 1: Covers positions around 1000
    a = pysam.AlignedSegment()
    a.query_name = "read_001"
    a.query_sequence = "A" * 50
    a.flag = 99
    a.reference_id = 0  # S288C_chrI
    a.reference_start = 980
    a.mapping_quality = 60
    a.cigar = ((0, 50),)  # 50M
    a.query_qualities = pysam.qualitystring_to_array("I" * 50)
    outf.write(a)

    # Read 2: Covers positions around 1500
    b = pysam.AlignedSegment()
    b.query_name = "read_002"
    b.query_sequence = "C" * 50
    b.flag = 147
    b.reference_id = 0
    b.reference_start = 1480
    b.mapping_quality = 59
    b.cigar = ((0, 50),)
    b.query_qualities = pysam.qualitystring_to_array("H" * 50)
    outf.write(b)

    # Read 3: Covers positions around 2000
    c = pysam.AlignedSegment()
    c.query_name = "read_003"
    c.query_sequence = "G" * 50
    c.flag = 83
    c.reference_id = 0
    c.reference_start = 1980
    c.mapping_quality = 58
    c.cigar = ((0, 50),)
    c.query_qualities = pysam.qualitystring_to_array("J" * 50)
    outf.write(c)

    # Read 4: Covers positions around 2958 (original test position)
    d = pysam.AlignedSegment()
    d.query_name = "read_004"
    d.query_sequence = "T" * 50
    d.flag = 99
    d.reference_id = 0
    d.reference_start = 2940
    d.mapping_quality = 60
    d.cigar = ((0, 50),)
    d.query_qualities = pysam.qualitystring_to_array("I" * 50)
    outf.write(d)

    # Read 5: Covers positions around 2965 (original test position)
    e = pysam.AlignedSegment()
    e.query_name = "read_005"
    e.query_sequence = "A" * 50
    e.flag = 147
    e.reference_id = 0
    e.reference_start = 2950
    e.mapping_quality = 60
    e.cigar = ((0, 50),)
    e.query_qualities = pysam.qualitystring_to_array("I" * 50)
    outf.write(e)

    # Read 6: On chrII at position 1000
    f = pysam.AlignedSegment()
    f.query_name = "read_006"
    f.query_sequence = "C" * 50
    f.flag = 99
    f.reference_id = 1  # S288C_chrII
    f.reference_start = 980
    f.mapping_quality = 60
    f.cigar = ((0, 50),)
    f.query_qualities = pysam.qualitystring_to_array("I" * 50)
    outf.write(f)

    # Read 7: On chrII at position 2000
    g = pysam.AlignedSegment()
    g.query_name = "read_007"
    g.query_sequence = "G" * 50
    g.flag = 147
    g.reference_id = 1
    g.reference_start = 1980
    g.mapping_quality = 60
    g.cigar = ((0, 50),)
    g.query_qualities = pysam.qualitystring_to_array("I" * 50)
    outf.write(g)

# Index the BAM
print("  Indexing BAM file...")
pysam.index(str(bam_path))

# ============================================================================
# 2. Create positions.txt file
# ============================================================================
print("  Creating positions.txt...")

pos_content = """# Example positions file
# Format: chrom pos (1-based)
S288C_chrI	1000
S288C_chrI	1500
S288C_chrI	2000
S288C_chrI	2958
S288C_chrI	2965
S288C_chrII	1000
S288C_chrII	2000
"""
pos_path.write_text(pos_content)

# ============================================================================
# 3. Create valid VCF file
# ============================================================================
print("  Creating example.vcf...")

vcf_content = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=S288C_chrI,length=10000>
##contig=<ID=S288C_chrII,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
S288C_chrI	1000	.	A	G	30	PASS	DP=50	GT	0/1
S288C_chrI	1500	.	C	T	35	PASS	DP=55	GT	0/1
S288C_chrI	2000	.	G	A	40	PASS	DP=60	GT	0/1
S288C_chrI	2958	.	A	G	30	PASS	DP=50	GT	0/1
S288C_chrI	2965	.	C	T	30	PASS	DP=50	GT	0/1
S288C_chrII	1000	.	T	C	30	PASS	DP=50	GT	0/1
S288C_chrII	2000	.	G	A	30	PASS	DP=50	GT	0/1
"""
vcf_path.write_text(vcf_content)

# ============================================================================
# 4. Create malformed VCF file for error testing
# ============================================================================
print("  Creating malformed.vcf...")

malformed_vcf_content = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##contig=<ID=S288C_chrI,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
S288C_chrI	1000	.	A	G	30	PASS	DP=50	GT	0/1
S288C_chrI	abc	.	C	T	30	PASS	DP=50	GT	0/1
S288C_chrI	1500	.	C	T	30	PASS	DP=50	GT	0/1
S288C_chrI	0	.	G	A	30	PASS	DP=50	GT	0/1
S288C_chrI	-100	.	T	C	30	PASS	DP=50	GT	0/1
S288C_chrI	2000	.	A	G	10	LowQual	DP=5	GT	0/1
S288C_chrI	2100	.	AT	CG	30	PASS	DP=50	GT	0/1
S288C_chrI	2200	.	A	G,T	30	PASS	DP=50	GT	0/1
S288C_chrI	2300	.	A		30	PASS	DP=50	GT	0/1
S288C_chrI
S288C_chrI	2958	.	A	G	30	PASS	DP=50	GT	0/1
S288C_chrI	2965	.	C	T	30	PASS	DP=50	GT	0/1
"""
malformed_vcf_path.write_text(malformed_vcf_content)

# ============================================================================
# Print summary
# ============================================================================
print("\n" + "="*70)
print("Example data created successfully in ./example_data/")
print("="*70)
print("\nFiles created:")
print(f"  1. {bam_path} (BAM file with 7 reads)")
print(f"  2. {bam_path}.bai (BAM index)")
print(f"  3. {pos_path} (7 positions in 2-column format)")
print(f"  4. {vcf_path} (7 valid SNP positions)")
print(f"  5. {malformed_vcf_path} (malformed VCF for error testing)")

print("\nTest commands:")
print("\n# Test with positions file:")
print(f"  python bamCol.py {bam_path} --pos-file {pos_path} --out results.csv")

print("\n# Test with VCF file:")
print(f"  python bamCol.py {bam_path} --vcf-file {vcf_path} --out results.csv")

print("\n# Test with malformed VCF (should show warnings):")
print(f"  python bamCol.py {bam_path} --vcf-file {malformed_vcf_path} --out results.csv")

print("\n# Test with region filtering:")
print(f"  python bamCol.py {bam_path} --vcf-file {vcf_path} \\")
print(f"    --exclude-region S288C_chrI:1400-1600 --out results.csv")

print("\n# Test with manual position:")
print(f"  python bamCol.py {bam_path} --pos S288C_chrI:2958 --out results.csv")

print("\n# Test with multiprocessing:")
print(f"  python bamCol.py {bam_path} --vcf-file {vcf_path} --process 2 --out results.csv")

print("\n" + "="*70)
