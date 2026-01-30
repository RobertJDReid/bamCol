FROM python:3.10-slim

# Metadata
LABEL maintainer="Robert J. D. Reid"
LABEL version="0.3.2"
LABEL description="bamCol.py - Extract per-read base calls from BAM files"

# Install system dependencies required for pysam
RUN apt-get update && apt-get install -y \
    gcc \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pysam
RUN pip install --no-cache-dir pysam

# Create app directory
WORKDIR /app

# Copy the script
COPY bamCol.py /app/

# Set the entrypoint to run bamCol.py
ENTRYPOINT ["python", "/app/bamCol.py"]

# Default command shows help if no arguments provided
CMD ["--help"]
