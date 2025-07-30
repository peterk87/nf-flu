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
PROFILE="test_assemblies,docker"
FASTAS_DIR=""
OUTDIR="nf-flu-results-assemblies"

print_help() {
    echo "Usage: $0 [-w WORKFLOW_PATH] [-m MEMORY] [-c CPU]"
    echo
    echo "Options:"
    echo "  -w WORKFLOW_PATH    Path to the Nextflow workflow (default: ${WORKFLOW_PATH})"
    echo "  -o OUTDIR           Path to nf-flu output directory (default: ${OUTDIR})"
    echo "  -f FASTAS_DIR       Path to directory with Influenza FASTA files (default: $FASTAS_DIR)"
    echo "  -m MEMORY           Memory allocation for the Nextflow run (default: ${MEMORY})"
    echo "  -c CPU              CPU allocation for the Nextflow run (default: ${CPU})"
    echo "  -p PROFILE          Nextflow profile to use (default: ${PROFILE})"
    echo "  -h                  Display this help message"
}

while getopts "w:o:m:c:p:h" opt; do
    case $opt in
        w) WORKFLOW_PATH=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        m) MEMORY=$OPTARG ;;
        c) CPU=$OPTARG ;;
        p) PROFILE=$OPTARG ;;
        f) FASTAS_DIR=$OPTARG ;;
        h) print_help; exit 0 ;;
        \?) error "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
        :) error "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
    esac
done

shift $((OPTIND-1))

if [[ "${1:-}" == "--" ]]; then
    shift
fi

info "Starting nf-flu assemblies test execution script with ${CPU} CPU cores and ${MEMORY} memory..."

VADR_MODEL_TARGZ_URL="https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flu/1.6.3-2/vadr-models-flu-1.6.3-2.tar.gz"
VADR_MODEL_TARGZ="vadr-models-flu-1.6.3-2.tar.gz"
FASTA_ZST_URL="https://api.figshare.com/v2/file/download/53449877"
CSV_ZST_URL="https://api.figshare.com/v2/file/download/53449874"
FASTA_ZST_FILE="influenza.fna.zst"
CSV_ZST_FILE="influenza.csv.zst"
FASTAS_URL="https://github.com/CFIA-NCFAD/nf-test-datasets/raw/refs/heads/nf-flu/fastas.tar.gz"

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

if [ -d "$FASTAS_DIR" ]; then
    info "Input directory with FASTA files specified: '$FASTAS_DIR'"
else
    info "Downloading default FASTA files from '$FASTAS_URL'"
    curl -SLk https://github.com/CFIA-NCFAD/nf-test-datasets/raw/refs/heads/nf-flu/fastas.tar.gz | tar -xzf -
    FASTAS_DIR="fastas/"
fi

info "Download FASTA and CSV files"
download_file "$FASTA_ZST_URL" "$FASTA_ZST_FILE"
download_file "$CSV_ZST_URL" "$CSV_ZST_FILE"

info "Download VADR model tar.gz"
download_file "$VADR_MODEL_TARGZ_URL" "$VADR_MODEL_TARGZ"

if [ -d "$WORKFLOW_PATH" ]; then
 info "Running Nextflow pipeline from local path: $WORKFLOW_PATH"
else
 info "Pulling Nextflow pipeline from remote: $WORKFLOW_PATH"
 nextflow pull "$WORKFLOW_PATH"
fi

nextflow run "$WORKFLOW_PATH" \
    -profile $PROFILE \
    -resume \
    --platform assemblies \
    --input $FASTAS_DIR \
    --outdir $OUTDIR \
    --ncbi_influenza_fasta $FASTA_ZST_FILE \
    --ncbi_influenza_metadata $CSV_ZST_FILE \
    --vadr_model_targz $VADR_MODEL_TARGZ \
    --max_cpus $CPU --max_memory "$MEMORY" $@
