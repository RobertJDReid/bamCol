# bam_col.py

**Extract per-read base calls at specific genomic positions from a BAM file using pysam.**

This script produces a CSV table listing the reads at the given positions in a reference genome, along with base calls, mapping quality, strand, and other information.
It is useful for inspecting allele composition, verifying variant calls, or extracting read-level evidence at specific coordinates.
The output file can be used to track SNP outcomes along specific reads providing evidence for strand exchanges.

---

## 🔧 Requirements

- Python 3.8+
- [pysam](https://pysam.readthedocs.io/en/latest/)

Installation:
```bash
mamba create -n bamcol python=3.10 pysam -y
mamba activate bamcol
```
or
```bash
pip install pysam
```

---

## 🚀 Usage

```bash
python bam_col.py <bamfile> [options]
```

Example:
```bash
python bam_col.py sample.bam --pos S288C_chrI:2941 --pos S288C_chrI:2947 --out alleles.csv
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

---

## 📄 Output Format

Each row in the output CSV represents a **read** overlapping a specified reference position.

| Column | Description |
|---------|--------------|
| `chrom` | Chromosome / reference name. |
| `pos` | 1-based reference position. |
| `read_id` | Read name (query name). |
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
python bam_col.py yeast.bam --pos S288C_chrIV:928598
```

**Include CIGAR strings:**
```bash
python bam_col.py yeast.bam --pos-file positions.txt --cigar --out calls_with_cigar.csv
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

---

## 🧠 Notes

- Requires an indexed BAM (`.bam.bai`).
- For large datasets, adjust `--max-depth` to control performance and memory use.

---

## 🪪 License

MIT License © 2025 Robert J. D. Reid  
Contributions welcome.
