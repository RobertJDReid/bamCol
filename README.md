# bamCol.py

_version 0.2.1_

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

Installation:

First, make a Python environment with pysam

```bash
conda create -n bamcol python=3.10 pysam -y
conda activate bamcol
```

Then `git clone repo`

subsitute the copied repository link for 'repo' in the command above

If pysam fails to install while creating an environment, you can activate the environment
and

```bash
pip install pysam
```

---

## 🐳 Docker Usage

### Building the Docker Image

```bash
docker build -t bam_col -f Dockerfile .
```

### Running with Docker

Run the container interactively with your data mounted:

```bash
docker run --rm -ti -v "$PWD:/app" bam_col
```

### Example Docker Command

Once inside the container:

```bash
python archive/bam_col.py example_data/example.bam  --vcf-file example_data/sk50_calls_nrp017_filtered.vcf.gz   --include-region-file example_data/regions.csv   --out alleles.csv 
```


---

## 📋 Required Inputs

### 1. BAM File (Required)
- **Format**: Binary Alignment Map (BAM)
- **Requirements**: 
  - Must be indexed (`.bam.bai` file in the same directory)
  - Create index with: `samtools index your_file.bam`
- **Example**: `sample.bam` with corresponding `sample.bam.bai`

### 2. Position Specification (At least one required)

You must provide positions using **at least one** of these methods:

#### Option A: Command-line positions (`--pos`)
- **Format**: `CHR:POS` (1-based coordinates)
- **Example**: `--pos S288C_chrI:2941 --pos S288C_chrI:2947`
- Can be repeated multiple times

#### Option B: Position file (`--pos-file`)
- **Format**: TSV or CSV file with columns: `chrom pos`
- **Requirements**:
  - Positions are 1-based
  - Header line must be commented with `#` or omitted
  - Lines starting with `#` are ignored
- **Example file** (`positions.csv`):
  ```
  #chrom,pos
  S288C_chrI,2941
  S288C_chrI,2947
  S288C_chrII,123927
  ```

#### Option C: VCF file (`--vcf-file`)
- **Format**: VCF/VCF.gz (Variant Call Format)
- **Requirements**:
  - Supports uncompressed `.vcf`, bgzip-compressed `.vcf.gz`, or regular gzip-compressed `.vcf.gz`
  - Only biallelic SNPs with FILTER=PASS (or empty filter) are used
  - VCF positions are automatically extracted (1-based)
- **Example**: `--vcf-file variants.vcf.gz`

#### Option D: Region filtering

Optionally combine with region filters to include/exclude specific genomic regions:

- **Include regions** (`--include-region` or `--include-region-file`):
  - Format: `CHR:START-END` (1-based, inclusive)
  - Example: `--include-region S288C_chrI:1000-5000`
  - Region file format: TSV/CSV with columns `chrom start end`

- **Exclude regions** (`--exclude-region`):
  - Format: `CHR:START-END` (1-based, inclusive)
  - Example: `--exclude-region S288C_chrI:2000-2100`

### 3. Output (Optional)
- **Default**: stdout (can be piped)
- **Specify file**: `--out output.csv`

---

## 🚀 Usage

```bash
python bamCol.py <bamfile> [options]
```

Examples:
```bash
python bamCol.py sample.bam --pos S288C_chrI:2941 --pos S288C_chrI:2947 --out alleles.csv
```

Works great with the multithreaded [pigz](https://zlib.net/pigz/) compression utility.

```bash
python bamCol.py sample.bam --pos-file SK1_SNPs.csv | pigz > zipped_output.csv.gz
```

Using VCF file with region filtering:
```bash
python bamCol.py sample.bam --vcf-file variants.vcf.gz \
  --include-region-file regions.csv --out alleles.csv
```

---

## ⚙️ Options

| Option | Description |
|--------|--------------|
| `bam` | Input BAM file (must be indexed `.bai`). |
| `--pos CHR:POS` | 1-based coordinate like `S288C_chrII:123927`. Can be repeated. |
| `--pos-file FILE` | TSV or CSV with `chrom pos` columns (1-based). Lines starting with `#` are ignored. |
| `--vcf-file FILE` | VCF/VCF.gz file. SNP positions will be extracted (supports bgzip or regular gzip). |
| `--include-region CHR:START-END` | Include only positions within region (repeatable, 1-based inclusive). |
| `--include-region-file FILE` | File with regions: `chrom start end` (1-based inclusive). |
| `--exclude-region CHR:START-END` | Exclude positions within region (repeatable, 1-based inclusive). |
| `--out FILE` | Output CSV file (default: stdout). |
| `--min-mapq N` | Minimum mapping quality (default: `0`). |
| `--min-bq N` | Minimum base quality (Phred, default: `0`). |
| `--max-depth N` | Maximum pileup depth (default: `100000`). |
| `--no-ignore-overlaps` | Count both mates if they overlap (default ignores overlaps to prevent double counting). |
| `--ignore-orphans` | Ignore reads whose mate is not properly paired. |
| `--cigar` | Include the read's CIGAR string in the output. Omit to exclude the column. |
| `--process n` | defines the number of system processes to use. Defaults to the total number of cores available from the OS. |
| `--include-secondary` | includes information from secondary alignments. |
| `--include-supplementary` | includes information from supplementary alignments. |

---

## 📄 Output Format

Each row in the output CSV represents a **read** overlapping a specified reference position.

| Column | Description |
|---------|--------------|
| `chrom` | Chromosome / reference name. |
| `pos` | 1-based reference position. |
| `read_id` | Read name (query name). |
| `read_pos` | Position on read aligned with current reference nucleotide |
| `call` | Base observed at that position (`A`, `T`, `G`, `C`, `DEL`, `REFSKIP`). |
| `is_del` | `True` if the read has a deletion at this position. |
| `is_refskip` | `True` if the read skips this reference position (e.g., spliced RNA-seq read). |
| `base_qual` | Phred base quality score. |
| `mapq` | Read mapping quality. |
| `strand` | `+` or `-` strand. |
| `cigar` | *(Optional)* CIGAR string (only if `--cigar` is used). |
| `flag` | SAM bitwise flag. |

---

## 🧩 Example

**Inspect base calls at a position:**
```bash
python bamCol.py yeast.bam --pos S288C_chrIV:928598
```

**Include CIGAR strings:**
```bash
python bamCol.py yeast.bam --pos-file positions.txt --cigar --out calls_with_cigar.csv
```

**Input position file (`positions.csv`):**
```
#chrom,pos
S288C_chrI,2941
S288C_chrI,2947
S288C_chrII,123927
```

note: `#` comments out the header line in the input file.
_The position file does not need a header and the script will fail unless commented out._

**Example output:**
```csv
chrom,pos,read_id,call,is_del,is_refskip,base_qual,mapq,strand,flag
S288C_chrI,2941,read_001,A,False,False,38,60,+,99
S288C_chrI,2941,read_002,DEL,True,False,,59,-,147
S288C_chrI,2947,read_003,G,False,False,41,60,+,83
```

With `--cigar`, an extra `cigar` column appears:
```csv
chrom,pos,read_id,call,is_del,is_refskip,base_qual,mapq,strand,cigar,flag
S288C_chrI,2941,read_001,A,False,False,38,60,+,76M,99
```

Installation can be tested by running:

```python make_example_data.py```

This script makes an example data folder
containing a very small bam file and positions file.
Then run:

```python bamCol.py example_data/example.bam --pos-file example_data/positions.csv```

Test output will be the following:

|chrom|pos|read_id|call|is_del|is_refskip|base_qual|mapq|strand|flag|
|---------|-----|--------|:--:|:--:|:--:|:--:|:--:|:--:|--:|
|S288C_chrI|2958|read_001|A|0|0|40|60|+|99|
|S288C_chrI|2958|read_002|C|0|0|39|59|-|147|
|S288C_chrI|2958|read_003|G|0|0|41|58|-|83|
|S288C_chrI|2965|read_001|A|0|0|40|60|+|99|
|S288C_chrI|2965|read_002|C|0|0|39|59|-|147|
|S288C_chrI|2965|read_003|G|0|0|41|58|-|83|

---

## 🧠 Notes

- Requires an indexed BAM (`.bam.bai`).
- For large datasets, adjust `--max-depth` to control performance and memory use.
- For large outputs, pipe into a compression utility such as `gzip` or `pigz`.

---

## 🪪 License

MIT License © 2025 Robert J. D. Reid  
Contributions welcome.
