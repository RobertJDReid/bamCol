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

---

## ⚙️ Options

| Option | Description |
|--------|--------------|
| `bam` | Input BAM file (must be indexed `.bai`). |
| `--pos CHR:POS` | 1-based coordinate like `S288C_chrII:123927`. Can be repeated. |
| `--pos-file FILE` | TSV or CSV with `chrom pos` columns (1-based). Lines starting with `#` are ignored. |
| `--out FILE` | Output CSV file (default: stdout). |
| `--min-mapq N` | Minimum mapping quality (default: `0`). |
| `--min-bq N` | Minimum base quality (Phred, default: `0`). |
| `--max-depth N` | Maximum pileup depth (default: `100000`). |
| `--no-ignore-overlaps` | Count both mates if they overlap (default ignores overlaps to prevent double counting). |
| `--ignore-orphans` | Ignore reads whose mate is not properly paired. |
| `--cigar` | Include the read’s CIGAR string in the output. Omit to exclude the column. |
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

### 🔎 VCF SNP Pre-filtering Helper

This repository includes an auxiliary R script for **pre-filtering SNPs from bcftools-generated VCF files** prior to downstream read-level analysis with `bamCol`. The script applies a set of stringent, reproducible filters to identify **high-confidence SNP positions**, optionally exporting both a filtered VCF and a simple CHROM/POS coordinate list.

Specifically, the script:

* retains **biallelic SNPs only** (excludes indels and multi-allelic sites),
* enforces a **minimum QUAL threshold**,
* filters by **site depth relative to the median INFO/DP**,
* selects sites with **AC = 2** (homozygous ALT in single-sample bcftools VCFs),
* optionally excludes specified chromosomes (e.g. mitochondrial DNA).

The resulting coordinate CSV can be used as input to `bamCol`, enabling efficient extraction of per-read allele information at robust, fixed SNP sites.

---

### 🔁 Typical Workflow to make a SNP position file

A common use case is to identify high-confidence SNP markers from a bcftools VCF and then extract per-read allele calls at those positions using `bamCol`:

```bash
# 0) Map diverged parent strain reads to the reference parent

minimap2 -ax asm5 ref_genome.fasta diverged_genome_reads.fastq.gz \
  | samtools view -b -q 50 \
  | samtools sort \
  > sample.bam

samtools index sample.bam

# 1) Call variants with bcftools (example)
bcftools mpileup -Ou -f reference_genome.fa sample.bam \
  | bcftools call -mv -Oz -o sample.vcf.gz
bcftools index sample.vcf.gz

# 2) Filter to high-confidence SNP positions
Rscript vcf_filt.R --vcf_file sample.vcf.gz --csv_only

# 3) Extract per-read base calls at filtered SNP sites
python bamCol.py sample.bam \
  --pos-file sample_filtered.csv \
  --out allele_calls.csv
```

Note that the `vcf_filt.R` script requires the following packages available on [CRAN](https://cran.r-project.org/)

- tidyverse
- vcfR
- optparse

---

## 🪪 License

MIT License © 2025 Robert J. D. Reid  
Contributions welcome.
