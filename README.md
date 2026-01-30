# bamCol.py

_version 0.3.2_

![sequenceRain](images/Align_rain2.png)

**Extract per-read base calls at specific genomic positions from a BAM file using pysam.**

Given a chromosomal position or positions, the script extracts base calls from a bam file of mapped sequence reads.
Output is a CSV table listing the read IDs, base calls, mapping quality, strand, and other information.
A positions list is automatically divided among available processor cores for parallel processing.
`bamCol` is useful for inspecting allele composition, verifying variant calls, or extracting read-level
evidence at specific coordinates.
The output file can be used to track SNP outcomes along specific reads providing evidence for strand exchanges
or other chromosomal alterations.
By default, only primary alignment reads are returned.
Options are available to include secondary and supplemental read information.

---

## 🔧 Requirements

- Python 3.8+
- [pysam](https://pysam.readthedocs.io/en/latest/)

---

## 📦 Installation

### Option 1: Conda Environment (Recommended)

First, make a Python environment with pysam:

```bash
conda create -n bamcol python=3.10 pysam -y
conda activate bamcol
```

Then clone the repository:

```bash
git clone <repository-url>
cd bamCol
```

If pysam fails to install while creating an environment, you can activate the environment and:

```bash
pip install pysam
```

### Option 2: Docker Container

Docker provides an isolated, reproducible environment without needing to manage Python dependencies.

#### Building the Docker Image

From the repository directory containing the `Dockerfile`:

```bash
docker build -t bamcol:0.3.2 .
```

#### Example Dockerfile

The repository includes a `Dockerfile`. For reference, it looks like this:

```dockerfile
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pysam
RUN pip install --no-cache-dir pysam

# Copy the script
WORKDIR /app
COPY bamCol.py /app/

# Set the entrypoint
ENTRYPOINT ["python", "/app/bamCol.py"]
```

#### Docker Tips

**Create an alias for easier usage:**
```bash
# Add to your ~/.bashrc or ~/.zshrc
alias bamcol='docker run --rm -v $(pwd):/data bamcol:0.3.2'

# Then use like a normal command
bamcol /data/sample.bam --pos chr1:1000 --out /data/results.csv
```

**Testing the Docker installation:**
```bash
# First, create example data (mount current directory)
docker run --rm -v $(pwd):/data -w /data python:3.10-slim \
    bash -c "pip install pysam && python /data/make_example_data.py"

# Then test bamCol
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/example_data/example.bam \
    --pos-file /data/example_data/positions.txt
```

---

## 🚀 Usage

### Direct Python Usage

```bash
python bamCol.py <bamfile> [options]
```

**Examples:**
```bash
python bamCol.py sample.bam --pos S288C_chrI:2941 --pos S288C_chrI:2947 --out alleles.csv
```

Works great with the multithreaded [pigz](https://zlib.net/pigz/) compression utility:

```bash
python bamCol.py sample.bam --pos-file SK1_SNPs.csv | pigz > zipped_output.csv.gz
```

### Docker Usage

When using Docker, you need to mount your data directory and use paths relative to the mount point.

**Basic usage with Docker:**
```bash
docker run --rm -v /path/to/data:/data bamcol:0.3.2 \
    /data/sample.bam --pos chr1:1000 --out /data/results.csv
```

**Using a position file:**
```bash
docker run --rm -v /path/to/data:/data bamcol:0.3.2 \
    /data/sample.bam --pos-file /data/positions.txt --out /data/results.csv
```

**Using a VCF file:**
```bash
docker run --rm -v /path/to/data:/data bamcol:0.3.2 \
    /data/sample.bam --vcf-file /data/variants.vcf.gz --out /data/results.csv
```

**With region filtering:**
```bash
docker run --rm -v /path/to/data:/data bamcol:0.3.2 \
    /data/sample.bam \
    --vcf-file /data/variants.vcf.gz \
    --exclude-region chr1:1000-2000 \
    --out /data/results.csv
```

**Using multiprocessing:**
```bash
docker run --rm -v /path/to/data:/data bamcol:0.3.2 \
    /data/sample.bam \
    --vcf-file /data/variants.vcf.gz \
    --process 4 \
    --out /data/results.csv
```

**Output to stdout (pipe to file or another tool):**
```bash
docker run --rm -v /path/to/data:/data bamcol:0.3.2 \
    /data/sample.bam --pos chr1:1000 > results.csv
```

**Docker volume mounting notes:**
- `-v /path/to/data:/data` mounts your local directory to `/data` in the container
- BAM files, position files, VCF files, and output files must all be in the mounted directory
- Use `/data/` prefix for all file paths inside the container

---

## ⚙️ Options

### Required
| Option | Description |
|--------|--------------|
| `bam` | Input BAM file (must be indexed `.bai`). |

### Position Input (at least one required)
| Option | Description |
|--------|--------------|
| `--pos CHR:POS` | 1-based coordinate like `S288C_chrII:123927`. Can be repeated. |
| `--pos-file FILE` | TSV or CSV with `chrom pos` columns (1-based). Lines starting with `#` are ignored. |
| `--vcf-file FILE` | VCF file. Biallelic SNP positions with FILTER=PASS will be extracted. |

### Region Filtering (optional)
| Option | Description |
|--------|--------------|
| `--include-region CHR:START-END` | Include only positions within region (repeatable). Format: `chr1:1000-2000` (1-based, inclusive). |
| `--exclude-region CHR:START-END` | Exclude positions within region (repeatable). Format: `chr1:1000-2000` (1-based, inclusive). |

### Output Options
| Option | Description |
|--------|--------------|
| `--out FILE` | Output CSV file (default: stdout). |

### Quality Filters
| Option | Description |
|--------|--------------|
| `--min-mapq N` | Minimum mapping quality (default: `0`). |
| `--min-bq N` | Minimum base quality (Phred, default: `0`). |
| `--max-depth N` | Maximum pileup depth (default: `100000`). |

### Alignment Options
| Option | Description |
|--------|--------------|
| `--no-ignore-overlaps` | Count both mates if they overlap (default ignores overlaps to prevent double counting). |
| `--ignore-orphans` | Ignore reads whose mate is not properly paired. |
| `--include-secondary` | Include information from secondary alignments. |
| `--include-supplementary` | Include information from supplementary alignments. |

### Performance & Output Format
| Option | Description |
|--------|--------------|
| `--process N` | Number of worker processes (default: all available CPU cores). |
| `--cigar` | Include the read's CIGAR string in the output. |
| `--version` | Show version number and exit. |

---

## 📄 Output Format

Each row in the output CSV represents a **read** overlapping a specified reference position.

| Column | Description |
|---------|--------------|
| `chrom` | Chromosome / reference name. |
| `pos` | 1-based reference position. |
| `read_id` | Read name (query name). |
| `read_pos` | Position on read aligned with current reference nucleotide. |
| `call` | Base observed at that position (`A`, `T`, `G`, `C`, `DEL`, `REFSKIP`). |
| `is_del` | `True` if the read has a deletion at this position. |
| `is_refskip` | `True` if the read skips this reference position (e.g., spliced RNA-seq read). |
| `base_qual` | Phred base quality score. |
| `mapq` | Read mapping quality. |
| `strand` | `+` or `-` strand. |
| `cigar` | *(Optional)* CIGAR string (only if `--cigar` is used). |
| `flag` | SAM bitwise flag. |

---

## 🧩 Examples

### Basic Usage

**Inspect base calls at a position:**
```bash
python bamCol.py yeast.bam --pos S288C_chrIV:928598
```

**Multiple positions:**
```bash
python bamCol.py sample.bam \
    --pos S288C_chrI:2941 \
    --pos S288C_chrI:2947 \
    --out alleles.csv
```

**Using a position file:**
```bash
python bamCol.py yeast.bam --pos-file positions.txt --out calls.csv
```

### VCF Input

**Extract positions from a VCF file:**
```bash
python bamCol.py sample.bam --vcf-file variants.vcf.gz --out allele_calls.csv
```

Note: Only biallelic SNPs with FILTER=PASS are used from VCF files.

### Region Filtering

**Exclude problematic regions:**
```bash
python bamCol.py sample.bam \
    --vcf-file variants.vcf.gz \
    --exclude-region chr1:1000-2000 \
    --exclude-region chr1:5000-6000 \
    --out filtered_calls.csv
```

**Focus on a specific region:**
```bash
python bamCol.py sample.bam \
    --vcf-file variants.vcf.gz \
    --include-region chr1:1000000-2000000 \
    --out targeted_calls.csv
```

**Combine include and exclude:**
```bash
python bamCol.py sample.bam \
    --pos-file positions.txt \
    --include-region chr1:1-10000000 \
    --exclude-region chr1:5000000-6000000 \
    --out calls.csv
```

### Advanced Options

**Include CIGAR strings:**
```bash
python bamCol.py yeast.bam --pos-file positions.txt --cigar --out calls_with_cigar.csv
```

**Use multiprocessing:**
```bash
python bamCol.py sample.bam --vcf-file variants.vcf.gz --process 8 --out calls.csv
```

**Apply quality filters:**
```bash
python bamCol.py sample.bam \
    --pos-file positions.txt \
    --min-mapq 30 \
    --min-bq 20 \
    --out high_quality_calls.csv
```

**Include secondary/supplementary alignments:**
```bash
python bamCol.py sample.bam \
    --pos-file positions.txt \
    --include-secondary \
    --include-supplementary \
    --out all_alignments.csv
```

---

## 📝 Input File Formats

### Position File Format

Tab or comma-separated file with columns: `chrom pos` (1-based).

**Example (`positions.txt`):**
```
#chrom,pos
S288C_chrI,2941
S288C_chrI,2947
S288C_chrII,123927
```

Note: Lines starting with `#` are treated as comments. The position file does not need a header.

### VCF File Format

Standard VCF format (version 4.x). Can be compressed (.vcf.gz) or uncompressed (.vcf).

**Filtering behavior:**
- Only biallelic SNPs are extracted (REF and ALT must be single bases)
- Only variants with FILTER=PASS (or missing filter)
- INDELs, multi-allelic sites, and non-PASS variants are skipped

**Example VCF:**
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000	.	A	G	30	PASS	DP=50
chr1	2000	.	C	T	35	PASS	DP=55
```

---

## 🧪 Testing Installation

Installation can be tested by running:

```bash
python make_example_data.py
```

This script creates an `example_data` folder containing:
- A small BAM file with 7 reads
- BAM index file
- A positions file
- A valid VCF file
- A malformed VCF file (for error testing)

Then test the installation:

```bash
python bamCol.py example_data/example.bam --pos-file example_data/positions.txt
```

**Expected output:**

|chrom|pos|read_id|read_pos|call|is_del|is_refskip|base_qual|mapq|strand|flag|
|---------|-----|--------|----|:--:|:--:|:--:|:--:|:--:|:--:|--:|
|S288C_chrI|1000|read_001|20|A|False|False|40|60|+|99|
|S288C_chrI|1500|read_002|20|C|False|False|39|59|-|147|
|S288C_chrI|2000|read_003|20|G|False|False|41|58|-|83|
|...|...|...|...|...|...|...|...|...|...|...|

Additional tests:

```bash
# Test VCF input
python bamCol.py example_data/example.bam --vcf-file example_data/example.vcf

# Test region filtering
python bamCol.py example_data/example.bam \
    --vcf-file example_data/example.vcf \
    --exclude-region S288C_chrI:1400-1600

# Test error handling with malformed VCF
python bamCol.py example_data/example.bam --vcf-file example_data/malformed.vcf
```

---

## 🧠 Notes

- Requires an indexed BAM (`.bam.bai`).
- For large datasets, adjust `--max-depth` to control performance and memory use.
- For large outputs, pipe into a compression utility such as `gzip` or `pigz`.
- Region filtering is optimized for 1-10 regions. For complex region filtering (>10 regions), consider pre-filtering your VCF with `bcftools` for better performance.
- When using `--vcf-file`, only high-quality biallelic SNPs are automatically extracted. For other variant types or failed-filter variants, use `--pos-file` instead.

---

## 🔄 Version History

### v0.3.2 (Current)
- Added minimal region filtering with `--include-region` and `--exclude-region`
- Performance warning for >10 regions
- Enhanced logging with region filtering statistics

### v0.3.1
- Improved error handling for malformed VCF and position files
- Added warnings for invalid positions and malformed data
- Better error messages throughout

### v0.3.0
- Added `--vcf-file` option for direct VCF input
- Automatic filtering of biallelic SNPs with FILTER=PASS
- Support for both bgzip and regular gzip compressed VCF files
- Position deduplication across all input sources

### v0.2.1
- Multiprocessing support with `--process` option
- Added `read_pos` column to output
- Optional CIGAR output with `--cigar` flag
- `--include-secondary` and `--include-supplementary` options
- Version tracking

---

## 📊 Typical Workflow

A common use case is to identify high-confidence SNP markers from a bcftools VCF and then extract per-read allele calls at those positions using `bamCol`:

```bash
# 0) Map diverged parent strain reads to the reference parent
minimap2 -ax asm5 ref_genome.fasta diverged_genome_reads.fastq.gz \
  | samtools view -b -q 50 \
  | samtools sort \
  > sample.bam

samtools index sample.bam

# 1) Call variants with bcftools
bcftools mpileup -Ou -f reference_genome.fa sample.bam \
  | bcftools call -mv -Oz -o sample.vcf.gz
bcftools index sample.vcf.gz

# 2) Extract per-read base calls at variant sites
python bamCol.py sample.bam \
  --vcf-file sample.vcf.gz \
  --out allele_calls.csv
```

### Alternative: Pre-filter VCF for complex scenarios

For complex filtering needs, use `bcftools` before bamCol:

```bash
# Filter VCF to specific regions using BED file
bcftools view -T regions.bed sample.vcf.gz -O z -o filtered.vcf.gz

# Extract per-read calls from filtered positions
python bamCol.py sample.bam --vcf-file filtered.vcf.gz --out allele_calls.csv
```

---

## 🔬 VCF SNP Pre-filtering Helper

This repository includes an auxiliary R script (`vcf_filt.R`) for **pre-filtering SNPs from bcftools-generated VCF files** prior to downstream read-level analysis with `bamCol`. The script applies a set of stringent, reproducible filters to identify **high-confidence SNP positions**, optionally exporting both a filtered VCF and a simple CHROM/POS coordinate list.

Specifically, the script:

* retains **biallelic SNPs only** (excludes indels and multi-allelic sites),
* enforces a **minimum QUAL threshold**,
* filters by **site depth relative to the median INFO/DP**,
* selects sites with **AC = 2** (homozygous ALT in single-sample bcftools VCFs),
* optionally excludes specified chromosomes (e.g. mitochondrial DNA).

The resulting coordinate CSV can be used as input to `bamCol` with `--pos-file`, enabling efficient extraction of per-read allele information at robust, fixed SNP sites.

**Required R packages** (available on [CRAN](https://cran.r-project.org/)):
- tidyverse
- vcfR
- optparse

**Note:** With the addition of `--vcf-file` support in bamCol v0.3.0+, many users can now skip this R pre-filtering step and use `--vcf-file` directly, as bamCol automatically filters for biallelic SNPs with FILTER=PASS. However, `vcf_filt.R` remains useful for more stringent custom filtering criteria.

---

## 🪪 License

MIT License © 2025 Robert J. D. Reid  
Contributions welcome.
