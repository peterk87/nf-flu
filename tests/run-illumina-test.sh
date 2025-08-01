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
SUBSAMPLE_INFLUENZA_FASTA=false
CPU=$(nproc)
MEMORY="8 GB"
PROFILE="test_illumina,docker"
OUTDIR="nf-flu-results-illumina"

print_help() {
    echo "Usage: $0 [-w WORKFLOW_PATH] [-m MEMORY] [-c CPU]"
    echo
    echo "Run the nf-flu test_illumina profile."
    echo
    echo "Options:"
    echo "  -w WORKFLOW_PATH              Path to the Nextflow workflow (default: ${WORKFLOW_PATH})"
    echo "  -o OUTDIR                     Path to nf-flu output directory (default: ${OUTDIR})"
    echo "  -s SUBSAMPLE_INFLUENZA_FASTA  Subsample the Influenza FASTA file the same way as GitHub Actions CI (default: ${SUBSAMPLE_INFLUENZA_FASTA})"
    echo "  -m MEMORY                     Memory allocation for the Nextflow run (default: ${MEMORY})"
    echo "  -c CPU                        CPU allocation for the Nextflow run (default: ${CPU})"
    echo "  -p PROFILE                    Nextflow profile to use (default: ${PROFILE})"
    echo "  -h                            Display this help message"
}

while getopts "w:o:s:m:c:p:h" opt; do
    case $opt in
        w) WORKFLOW_PATH=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        s) SUBSAMPLE_INFLUENZA_FASTA=true ;;
        m) MEMORY=$OPTARG ;;
        c) CPU=$OPTARG ;;
        p) PROFILE=$OPTARG ;;
        h) print_help; exit 0 ;;
        \?) error "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
        :) error "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
    esac
done

shift $((OPTIND-1))

if [[ "${1:-}" == "--" ]]; then
    shift
fi

info "Starting nf-flu Illumina test execution script with ${CPU} CPU cores and ${MEMORY} memory..."

VADR_MODEL_TARGZ_URL="https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flu/1.6.3-2/vadr-models-flu-1.6.3-2.tar.gz"
VADR_MODEL_TARGZ="vadr-models-flu-1.6.3-2.tar.gz"
FASTA_ZST_URL="https://api.figshare.com/v2/file/download/53449877"
CSV_ZST_URL="https://api.figshare.com/v2/file/download/53449874"
FASTA_ZST_FILE="influenza.fna.zst"
CSV_ZST_FILE="influenza.csv.zst"

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

info "Download FASTA and CSV files"
download_file "$FASTA_ZST_URL" "$FASTA_ZST_FILE"
download_file "$CSV_ZST_URL" "$CSV_ZST_FILE"

if [[ $SUBSAMPLE_INFLUENZA_FASTA ]] ; then
 info "Subsampling Influenza FASTA with 'seqtk sample -s 789'"
 # curl --silent -SLk ${FASTA_ZST_URL} | zstdcat | seqtk sample -s 789 - 10000 | zstd -ck > influenza-10k.fna.zst
 SEQTK_DOCKER_IMAGE="quay.io/biocontainers/seqtk:1.4--he4a0461_2"
 info "Pulling seqtk Docker image from $SEQTK_DOCKER_IMAGE"
 docker pull $SEQTK_DOCKER_IMAGE
 zstdcat ${FASTA_ZST_FILE} | docker run -i --rm $SEQTK_DOCKER_IMAGE seqtk sample -s 789 - 10000 | zstd -ck > influenza-10k.fna.zst
 FASTA_ZST_FILE="influenza-10k.fna.zst"
fi

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
    --outdir $OUTDIR \
    --ncbi_influenza_fasta $FASTA_ZST_FILE \
    --ncbi_influenza_metadata $CSV_ZST_FILE \
    --vadr_model_targz $VADR_MODEL_TARGZ \
    --max_cpus $CPU --max_memory "$MEMORY" $@
