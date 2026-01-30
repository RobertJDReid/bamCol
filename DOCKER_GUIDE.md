# Docker Quick Reference for bamCol.py

## Building the Image

```bash
# Build from repository directory
docker build -t bamcol:0.3.2 .

# Tag as latest
docker tag bamcol:0.3.2 bamcol:latest

# Verify the image
docker images | grep bamcol
```

---

## Basic Usage

### Simple position query
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --pos chr1:1000 --out /data/results.csv
```

### Using a position file
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --pos-file /data/positions.txt --out /data/results.csv
```

### Using a VCF file
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --vcf-file /data/variants.vcf.gz --out /data/results.csv
```

---

## Advanced Usage

### With region filtering
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam \
    --vcf-file /data/variants.vcf.gz \
    --exclude-region chr1:1000-2000 \
    --exclude-region chr2:5000-6000 \
    --out /data/results.csv
```

### With multiprocessing
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam \
    --vcf-file /data/variants.vcf.gz \
    --process 8 \
    --out /data/results.csv
```

### With quality filters
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam \
    --pos-file /data/positions.txt \
    --min-mapq 30 \
    --min-bq 20 \
    --out /data/results.csv
```

### Output to stdout (pipe to file)
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --pos chr1:1000 > results.csv
```

### Compress output on the fly
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --pos-file /data/positions.txt | gzip > results.csv.gz
```

---

## Volume Mounting

### Mount current directory
```bash
-v $(pwd):/data
```
Use `/data/` prefix for all files in the container.

### Mount specific directory
```bash
-v /path/to/bam/files:/bam -v /path/to/output:/output
```
Then use:
- `/bam/sample.bam` for input
- `/output/results.csv` for output

### Mount multiple directories
```bash
docker run --rm \
    -v /data/bams:/bams \
    -v /data/vcfs:/vcfs \
    -v /data/results:/results \
    bamcol:0.3.2 \
    /bams/sample.bam \
    --vcf-file /vcfs/variants.vcf.gz \
    --out /results/calls.csv
```

---

## Creating an Alias

Add to `~/.bashrc` or `~/.zshrc`:

```bash
alias bamcol='docker run --rm -v $(pwd):/data bamcol:0.3.2'
```

Then use like a regular command:
```bash
bamcol /data/sample.bam --pos chr1:1000 --out /data/results.csv
```

---

## Testing the Docker Image

### Create test data
```bash
# Download or create a small BAM file for testing
# Run the example data generator (if you have it mounted)
docker run --rm -v $(pwd):/data python:3.10-slim bash -c \
    "pip install pysam && python /data/make_example_data.py"
```

### Run test
```bash
docker run --rm -v $(pwd):/data bamcol:0.3.2 \
    /data/example_data/example.bam \
    --pos-file /data/example_data/positions.txt
```

### Check version
```bash
docker run --rm bamcol:0.3.2 --version
```

### View help
```bash
docker run --rm bamcol:0.3.2 --help
```

---

## Troubleshooting

### Permission issues
If you get permission errors with output files:
```bash
# Run with your user ID
docker run --rm --user $(id -u):$(id -g) \
    -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --pos chr1:1000 --out /data/results.csv
```

### File not found errors
Make sure:
1. The file exists in the mounted directory
2. You're using the correct path inside the container (e.g., `/data/file.bam`)
3. The BAM index file (`.bai`) is in the same directory as the BAM

### Container won't start
```bash
# Check if image exists
docker images | grep bamcol

# Rebuild if necessary
docker build -t bamcol:0.3.2 .

# Check Docker is running
docker ps
```

---

## Advanced Docker Options

### Limit CPU usage
```bash
docker run --rm --cpus=4 -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --vcf-file /data/variants.vcf.gz --process 4
```

### Limit memory
```bash
docker run --rm --memory=8g -v $(pwd):/data bamcol:0.3.2 \
    /data/sample.bam --vcf-file /data/variants.vcf.gz
```

### Run interactively (for debugging)
```bash
docker run -it --rm -v $(pwd):/data --entrypoint bash bamcol:0.3.2
# Inside container:
python /app/bamCol.py /data/sample.bam --pos chr1:1000
```

---

## Pushing to Docker Hub (Optional)

If you want to share the image:

```bash
# Tag with your Docker Hub username
docker tag bamcol:0.3.2 yourusername/bamcol:0.3.2
docker tag bamcol:0.3.2 yourusername/bamcol:latest

# Login to Docker Hub
docker login

# Push
docker push yourusername/bamcol:0.3.2
docker push yourusername/bamcol:latest
```

Then others can use it:
```bash
docker pull yourusername/bamcol:0.3.2
docker run --rm -v $(pwd):/data yourusername/bamcol:0.3.2 ...
```

---

## Multi-platform Build (ARM64 + AMD64)

For compatibility with both Intel/AMD and Apple Silicon:

```bash
# Create buildx builder
docker buildx create --name multiarch --use

# Build for multiple platforms
docker buildx build --platform linux/amd64,linux/arm64 \
    -t yourusername/bamcol:0.3.2 \
    --push .
```

---

## Clean Up

### Remove old containers
```bash
docker container prune
```

### Remove old images
```bash
docker image prune
```

### Remove specific image
```bash
docker rmi bamcol:0.3.2
```

### Remove all bamcol images
```bash
docker rmi $(docker images bamcol -q)
```
