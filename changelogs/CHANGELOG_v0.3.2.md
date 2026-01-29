# bamCol.py v0.3.2 - Minimal Region Filtering

## Changes from v0.3.1 to v0.3.2

### New Features

#### 1. **Region-Based Position Filtering**
Added minimal region filtering to handle common use cases without external tools.

**New command-line arguments:**

```bash
--include-region CHR:START-END
    Include only positions within the specified region (repeatable).
    Format: CHR:START-END (1-based, inclusive)
    Example: --include-region chr1:1000-2000
    
--exclude-region CHR:START-END
    Exclude positions within the specified region (repeatable).
    Format: CHR:START-END (1-based, inclusive)
    Example: --exclude-region chr1:5000-6000
```

**Key characteristics:**
- Both arguments are repeatable (can specify multiple regions)
- Supports both `-` and `..` as separators (chr1:1000-2000 or chr1:1000..2000)
- 1-based coordinates (inclusive on both ends)
- Auto-swaps start/end if provided in wrong order
- Validates that positions are >= 1
- Recommended for up to ~10 regions (see performance notes)

---

### New Functions

#### `parse_region_token(token: str) -> Tuple[str, int, int]`
Parses region strings like "chr1:1000-2000" into (chromosome, start, end) tuples.

**Features:**
- Accepts both `-` and `..` as separators
- Auto-swaps if start > end
- Validates positions >= 1
- Clear error messages for invalid formats

**Examples:**
```python
parse_region_token("chr1:1000-2000")  # → ("chr1", 1000, 2000)
parse_region_token("chr1:2000-1000")  # → ("chr1", 1000, 2000) - auto-swapped
parse_region_token("chr1:1000..2000") # → ("chr1", 1000, 2000) - .. separator
```

#### `in_regions(chrom: str, pos: int, regions: List[Tuple[str, int, int]]) -> bool`
Checks if a position falls within any of the specified regions.

**Returns:** True if position is in any region, False otherwise

---

### Enhanced Functions

#### `read_positions()` - Updated signature
```python
def read_positions(
    pos_args: List[str],
    pos_file: Optional[str],
    vcf_file: Optional[str],
    include_regions: List[Tuple[str, int, int]],  # NEW
    exclude_regions: List[Tuple[str, int, int]]   # NEW
) -> List[Tuple[str, int]]:
```

**New behavior:**
1. Collects positions from all sources (--pos, --pos-file, --vcf-file)
2. Applies include region filter (if specified)
3. Applies exclude region filter (if specified)
4. Deduplicates and sorts

**Filtering logic:**
- If `include_regions` specified: position must be in at least one include region
- If `exclude_regions` specified: position must not be in any exclude region
- Both can be combined for complex filtering

---

### Usage Examples

#### Example 1: Exclude a problematic region
```bash
# Exclude positions in low-coverage region
python bamCol_v0.3.2.py sample.bam \
    --vcf-file variants.vcf.gz \
    --exclude-region chr1:1000-2000 \
    --out results.csv
```

#### Example 2: Focus on specific region
```bash
# Only analyze positions in target region
python bamCol_v0.3.2.py sample.bam \
    --vcf-file variants.vcf.gz \
    --include-region chr1:1000000-2000000 \
    --out results.csv
```

#### Example 3: Multiple exclusions
```bash
# Exclude multiple known bad regions
python bamCol_v0.3.2.py sample.bam \
    --pos-file positions.txt \
    --exclude-region chr1:1000-2000 \
    --exclude-region chr1:5000-6000 \
    --exclude-region chr2:10000-11000 \
    --out results.csv
```

#### Example 4: Combine include and exclude
```bash
# Analyze chromosome 1, but exclude centromere
python bamCol_v0.3.2.py sample.bam \
    --vcf-file variants.vcf.gz \
    --include-region chr1:1-150000000 \
    --exclude-region chr1:60000000-65000000 \
    --out results.csv
```

#### Example 5: Mix position sources with regions
```bash
# Combine manual positions, file, and VCF; exclude bad regions
python bamCol_v0.3.2.py sample.bam \
    --pos chr1:1000 \
    --pos chr1:2000 \
    --pos-file additional.txt \
    --vcf-file variants.vcf.gz \
    --exclude-region chr1:1500-1600 \
    --exclude-region chr2:5000-6000 \
    --process 4 \
    --out results.csv
```

---

### Enhanced Logging

New stderr output includes region filtering information:
```
Using 8 worker process(es) out of 8 available CPU(s).
Region filtering: 1 include, 2 exclude
Processing 1247 position(s).
Wrote 45821 records
```

---

### Performance Considerations

**Algorithm:** O(n × m) where n = positions, m = regions

**Benchmarks:**
- 1,000 positions × 5 regions = 5,000 comparisons ✅ Fast
- 10,000 positions × 10 regions = 100,000 comparisons ✅ Acceptable
- 100,000 positions × 100 regions = 10M comparisons ⚠️ Slow

**Recommendations:**
1. **For ≤10 regions:** Use `--include-region` and `--exclude-region` (optimal)
2. **For >10 regions:** Pre-filter with bcftools (much faster)
3. **For >100 regions:** Definitely use bcftools with indexed BED files

**Performance warning:**
If you specify >10 total regions, you'll see:
```
Warning: Using 15 regions may impact performance. 
For large region lists, consider pre-filtering with bcftools.
```

---

### Alternative: Using bcftools for Complex Filtering

For complex region filtering, bcftools is significantly faster:

```bash
# Create BED file with many regions
cat > exclude_regions.bed << EOF
chr1	1000	2000
chr1	5000	6000
chr2	10000	11000
... (hundreds more)
EOF

# Pre-filter VCF with bcftools (very fast with indexed BED)
bcftools view -T ^exclude_regions.bed variants.vcf.gz -O z -o filtered.vcf.gz
tabix -p vcf filtered.vcf.gz

# Then use bamCol with pre-filtered VCF
python bamCol_v0.3.2.py sample.bam --vcf-file filtered.vcf.gz --out results.csv
```

This approach is orders of magnitude faster for large region lists.

---

### Error Handling Improvements

**Invalid region format:**
```bash
$ python bamCol_v0.3.2.py test.bam --exclude-region "chr1:abc-2000"
error: Bad region 'chr1:abc-2000'. Expected format: CHR:START-END (e.g., chr1:1000-2000)
```

**Invalid positions in region:**
```bash
$ python bamCol_v0.3.2.py test.bam --exclude-region "chr1:0-2000"
error: Region positions must be >= 1: chr1:0-2000
```

**Auto-correction (no error):**
```bash
# Start > End: automatically swapped
--exclude-region chr1:2000-1000  # Becomes chr1:1000-2000
```

---

### Implementation Details

**Code additions:**
- `parse_region_token()`: 25 lines
- `in_regions()`: 8 lines
- `read_positions()`: Updated to accept and apply region filters
- `main()`: Region parsing and validation logic
- Total new code: ~60 lines

**Dependencies:**
- Added `import re` for region parsing regex
- No new external dependencies

---

### Testing Region Filtering

#### Test 1: Basic exclusion
```bash
# positions.csv contains: chr1:1000, chr1:1500, chr1:2000
python bamCol_v0.3.2.py test.bam \
    --pos-file positions.csv \
    --exclude-region chr1:1400-1600

# Expected: Processes chr1:1000 and chr1:2000 only (1500 excluded)
# Output: Processing 2 position(s).
```

#### Test 2: Basic inclusion
```bash
# positions.csv contains: chr1:1000, chr1:2000, chr2:3000
python bamCol_v0.3.2.py test.bam \
    --pos-file positions.csv \
    --include-region chr1:1-10000

# Expected: Processes chr1:1000 and chr1:2000 only (chr2:3000 excluded)
# Output: Processing 2 position(s).
```

#### Test 3: Combined filtering
```bash
python bamCol_v0.3.2.py test.bam \
    --vcf-file test.vcf \
    --include-region chr1:1-1000000 \
    --exclude-region chr1:50000-60000

# Expected: Only chr1 positions, excluding 50000-60000 range
```

---

## Migration Guide

### From v0.3.1 to v0.3.2

**No Breaking Changes** - All existing commands work identically.

**New capabilities:**
```bash
# Old way (external tool required):
bcftools view -T ^exclude.bed variants.vcf.gz > filtered.vcf
python bamCol_v0.3.1.py test.bam --vcf-file filtered.vcf

# New way (built-in for simple cases):
python bamCol_v0.3.2.py test.bam \
    --vcf-file variants.vcf.gz \
    --exclude-region chr1:1000-2000
```

---

## When to Use Region Filtering

### Use `--include-region` / `--exclude-region` when:
- ✅ You have 1-10 regions to filter
- ✅ Quick ad-hoc analysis
- ✅ Excluding known problematic regions
- ✅ Filtering non-VCF position sources
- ✅ Iterative analysis (tweaking regions)

### Use bcftools pre-filtering when:
- ✅ You have >10 regions
- ✅ Using indexed BED files
- ✅ Complex region logic
- ✅ Maximum performance needed
- ✅ Reusing filtered VCF multiple times

---

## Updated example data for testing installation

## Upgrade Recommendation

**Recommended for all users.**

v0.3.2 adds convenient region filtering for common workflows without adding complexity for users who don't need it. The feature is entirely optional and doesn't affect existing functionality.

**Best for:**
- Users who need to exclude 1-5 problematic genomic regions
- Quick targeted analysis without creating intermediate files
- Filtering positions from multiple sources (not just VCF)

---

## Version History

- **v0.3.2** (Current) - Added minimal region filtering
- **v0.3.1** - Bug fixes and improved error handling
- **v0.3.0** - Added VCF file support
- **v0.2.1** - Multiprocessing and enhanced features
