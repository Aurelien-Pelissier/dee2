#!/bin/bash

ORG="kphaffii"
TAX_ID="460519"

BASE_DIR="/home/user/x-alife/user/dee2"
LIST_FILE="${BASE_DIR}/${ORG}_accessions.txt"
RESULTS_DIR="${BASE_DIR}/results"
LOGS_DIR="${BASE_DIR}/logs"

# OPTIMIZED FOR n1-standard-96 (12 jobs x 8 threads fills all 96 cores perfectly)
MAX_JOBS=12

mkdir -p "$RESULTS_DIR"
mkdir -p "$LOGS_DIR"
mkdir -p "${BASE_DIR}/ref"

# 1. Download the list of SRA runs (if missing or empty)
if [ ! -f "$LIST_FILE" ] || [ ! -s "$LIST_FILE" ]; then
    echo "Fetching transcriptomic SRA accessions for $ORG from EBI ENA API..."
    curl -s "https://www.ebi.ac.uk/ena/portal/api/search?query=tax_tree(${TAX_ID})&result=read_run&fields=run_accession,library_source&format=tsv&limit=0" \
        | awk -F'\t' '$2 ~ /TRANSCRIPTOMIC/ {print $1}' > "$LIST_FILE"
    echo "Found $(wc -l < "$LIST_FILE") RNA-Seq runs to process."
fi

# 2. Define the core processing worker function
process_srr() {
    local SRR="$1"
    local ORG="$2"
    local RESULTS_DIR="$3"
    local BASE_DIR="$4"
    local LOGS_DIR="$5"

    # Define a clean target path to check for completed runs
    # This checks for either a custom named zip, a standard zip, or an unzipped result directory
    if [ -s "$RESULTS_DIR/${SRR}.${ORG}.zip" ] || [ -s "$RESULTS_DIR/${SRR}/${SRR}.zip" ] || [ -s "$RESULTS_DIR/${SRR}.zip" ]; then
        echo "Skipping $SRR (already successfully processed)."
        return 0
    fi

    # Staggered start up to 30 seconds to protect network I/O from getting throttled
    local DELAY=$(( (RANDOM % 30) + 1 ))
    sleep $DELAY

    echo "=== Starting background job for $SRR (stagger delay: ${DELAY}s) ==="

    # Create a dedicated, isolated output directory for this specific run on the host 
    # to prevent parallel workers from clobbering each other's files.
    local RUN_OUTPUT_DIR="${RESULTS_DIR}/${SRR}"
    mkdir -p "$RUN_OUTPUT_DIR"

    # Route internal container outputs directly to the host machine using a volume mount.
    # We patch the volunteer_pipeline.sh execution by ensuring the container outputs drop into /dee2/mnt
    docker run --name "dee2_${SRR}" \
        -v "${BASE_DIR}/pipeline/volunteer_pipeline.sh:/dee2/code/volunteer_pipeline.sh" \
        -v "${BASE_DIR}/pipeline/FastqPairer.pl:/dee2/code/FastqPairer.pl" \
        -v "${BASE_DIR}/ref:/dee2/ref" \
        -v "${RUN_OUTPUT_DIR}:/dee2/mnt" \
        mziemann/tallyup -s "$ORG" -a "$SRR" > "${LOGS_DIR}/${SRR}.log" 2>&1

    # Fix ownership of the output directory and its contents before host-level move/rm
    # This ensures the host user has permissions to move the zip and delete the temp folder.
    docker run --rm -v "${RUN_OUTPUT_DIR}:/mnt" alpine chown -R $(id -u):$(id -g) /mnt

    # Verify if the container successfully generated data into our host-mounted directory
    # (Checking if the directory exists and contains files)
    if [ -d "$RUN_OUTPUT_DIR" ] && [ "$(ls -A "$RUN_OUTPUT_DIR" 2>/dev/null)" ]; then
        echo "Successfully saved $SRR"
        
        # If the pipeline left a zip inside the folder, pull it up one level for clean aesthetics
        # volunteer_pipeline.sh typically creates ${SRR}.${ORG}.zip in /dee2/mnt
        if [ -f "$RUN_OUTPUT_DIR/${SRR}.${ORG}.zip" ]; then
            mv "$RUN_OUTPUT_DIR/${SRR}.${ORG}.zip" "$RESULTS_DIR/${SRR}.${ORG}.zip"
            # Clean up the now redundant subfolder if everything was inside that zip
            rm -rf "$RUN_OUTPUT_DIR"
        elif [ -f "$RUN_OUTPUT_DIR/${SRR}.zip" ]; then
            mv "$RUN_OUTPUT_DIR/${SRR}.zip" "$RESULTS_DIR/${SRR}.${ORG}.zip"
            rm -rf "$RUN_OUTPUT_DIR"
        fi

        # Delete successful logs to save storage space, keeping failed ones for troubleshooting
        rm -f "${LOGS_DIR}/${SRR}.log"
    else
        echo "Error processing $SRR (Check logs/${SRR}.log for details)."
        # Remove the empty run folder so it attempts a retry on the next pipeline execution
        rm -rf "$RUN_OUTPUT_DIR"
    fi

    # Cleanup container
    docker rm -v "dee2_${SRR}" >/dev/null 2>&1
}

# Export environment variables so the sub-shells spawned by parallel can access them
export -f process_srr
export ORG RESULTS_DIR BASE_DIR LOGS_DIR

echo "------------------------------------------------"
echo "Launching parallel bulk processing queue..."
echo "Concurrency limit: $MAX_JOBS workers | Total runs: $(wc -l < "$LIST_FILE")"
echo "------------------------------------------------"

# Wipe out any dead or lingering 'dee2_' containers to prevent name collisions
echo "Cleaning up any old, lingering containers..."
docker rm -f $(docker ps -a -q --filter name=dee2_) >/dev/null 2>&1

# 3. Stream the accessions into GNU Parallel to maximize CPU utilization smoothly
cat "$LIST_FILE" | parallel --ungroup -j "$MAX_JOBS" "bash -c 'process_srr \"{}\" \"$ORG\" \"$RESULTS_DIR\" \"$BASE_DIR\" \"$LOGS_DIR\"'"

echo "Bulk processing complete! Parallel execution queue has finished processing."
