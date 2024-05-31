#!/bin/bash

set -euo pipefail

error_handler() {
    echo -e "\n\033[1;31mError on line $1\033[0m"
    # Perform any cleanup or logging here
}

cleanup() {
    echo "Cleaning up before exiting..."
    # Cleanup commands go here
}

handle_interrupt() {
    echo -e "\n\033[1;31mERROR:\033[1m Script interrupted...\033[0m"
    cleanup
    exit 1
}

# Trap ERR signal to handle errors
trap 'error_handler $LINENO' ERR

# Trap EXIT signal to perform cleanup
trap cleanup EXIT

# Trap SIGINT and SIGTERM to handle interruptions
trap handle_interrupt SIGINT SIGTERM


error() {
  echo -e "$(date -Is)  \033[1;31mERROR: \033[0m\033[1m$1\033[0m"
}

info() {
  echo -e "$(date -Is)  \033[1;32mINFO: \033[0m\033[1m$1\033[0m"
}


WORKFLOW_PATH="CFIA-NCFAD/nf-flu"
CPU=$(nproc)
MEMORY="8 GB"

print_help() {
    echo "Usage: $0 [-w WORKFLOW_PATH] [-m MEMORY] [-c CPU]"
    echo
    echo "Options:"
    echo "  -w WORKFLOW_PATH    Path to the Nextflow workflow (default: ${WORKFLOW_PATH})"
    echo "  -m MEMORY           Memory allocation for the Nextflow run (default: ${MEMORY})"
    echo "  -c CPU              CPU allocation for the Nextflow run (default: ${CPU})"
    echo "  -h                  Display this help message"
}

while getopts "w:m:c:h" opt; do
    case $opt in
        w) WORKFLOW_PATH=$OPTARG ;;
        m) MEMORY=$OPTARG ;;
        c) CPU=$OPTARG ;;
        h) print_help; exit 0 ;;
        \?) error "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
        :) error "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
    esac
done

info "Starting nf-flu Nanopore test execution script with ${CPU} CPU cores and ${MEMORY} memory..."

FASTA_ZST_URL="https://api.figshare.com/v2/file/download/41415330"
CSV_ZST_URL="https://api.figshare.com/v2/file/download/41415333"

download_file() {
    local url=$1
    local output=$2

    if [ ! -f "$output" ]; then
        info "Downloading $output from $url..."
        curl --silent -SLk "$url" -o "$output"
        info "Downloaded $output."
    else
        info "$output already exists. Skipping download."
    fi
}

create_samplesheet() {
    local samplesheet=$1
    echo "sample,reads" | tee "$samplesheet"
    echo "ERR6359501-10k,$(realpath reads/ERR6359501-10k.fastq.gz)" | tee -a "$samplesheet"
    echo "ERR6359501,$(realpath reads/run1)" | tee -a "$samplesheet"
    echo "ERR6359501,$(realpath reads/run2)" | tee -a "$samplesheet"
    echo "SRR24826962,$(realpath reads/SRR24826962.fastq.gz)" | tee -a "$samplesheet"
    echo "ntc-bc15,$(realpath reads/ntc-bc15.fastq.gz)" | tee -a "$samplesheet"
    echo "ntc-bc31,$(realpath reads/ntc-bc31.fastq.gz)" | tee -a "$samplesheet"
    echo "ntc-bc47,$(realpath reads/ntc-bc47.fastq.gz)" | tee -a "$samplesheet"
}

FASTA_ZST_FILE="influenza.fna.zst"
CSV_ZST_FILE="influenza.csv.zst"

# Create directories
mkdir -p reads/{run1,run2}

# Download test data files
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/SRR24826962.sampled.fastq.gz" "reads/SRR24826962.fastq.gz"
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/ERR6359501-10k.fastq.gz" "reads/ERR6359501-10k.fastq.gz"
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/run1-s11-ERR6359501.fastq.gz" "reads/run1/s11-ERR6359501.fastq.gz"
if [ ! -f "reads/run1/s1-ERR6359501.fastq" ]; then
    download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/run1-s1-ERR6359501.fastq.gz" "reads/run1/s1-ERR6359501.fastq.gz"
    gunzip reads/run1/s1-ERR6359501.fastq.gz
else
    info "'reads/run1/s1-ERR6359501.fastq' already exists!"
fi
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/run2-s22-ERR6359501.fastq.gz" "reads/run2/s22-ERR6359501.fastq.gz"
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/run2-s2-ERR6359501.fastq.gz" "reads/run2/s2-ERR6359501.fastq.gz"
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/ntc-bc15.fastq.gz" "reads/ntc-bc15.fastq.gz"
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/ntc-bc31.fastq.gz" "reads/ntc-bc31.fastq.gz"
download_file "https://github.com/CFIA-NCFAD/nf-test-datasets/raw/nf-flu/nanopore/fastq/ntc-bc47.fastq.gz" "reads/ntc-bc47.fastq.gz"

info "Creating samplesheet.csv"
create_samplesheet "samplesheet.csv"

info "Download FASTA and CSV files"
download_file "$FASTA_ZST_URL" "$FASTA_ZST_FILE"
download_file "$CSV_ZST_URL" "$CSV_ZST_FILE"

if [ -d "$WORKFLOW_PATH" ]; then
 info "Running Nextflow pipeline from local path: $WORKFLOW_PATH"
else
 info "Pulling Nextflow pipeline from remote: $WORKFLOW_PATH"
 nextflow pull "$WORKFLOW_PATH"
fi

nextflow run "$WORKFLOW_PATH" \
    -profile test_nanopore,docker \
    -resume \
    --platform nanopore \
    --input samplesheet.csv \
    --ncbi_influenza_fasta $FASTA_ZST_FILE \
    --ncbi_influenza_metadata $CSV_ZST_FILE \
    --max_cpus $CPU --max_memory "$MEMORY"
