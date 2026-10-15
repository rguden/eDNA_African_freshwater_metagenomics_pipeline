#!/bin/bash -l

# ===============================================================================
# eDNA Illumina workflow: metagenomics pipeline from samples from African freshwater systems
# Author: Rodgee Mae Guden
# ===============================================================================

set -euo pipefail

# ----------------- CONFIG -----------------
READS_DIR="${READS_DIR:-/staging/leuven/stg_00100/Rodgee/caves_data/concatenated_raw_reads}"
METADATA="${METADATA:-/staging/leuven/stg_00100/Rodgee/caves_data/metadata/biological_metadata.tsv}"
OUTPUT_DIR="${OUTPUT_DIR:-/staging/leuven/stg_00100/Rodgee/caves_data/output}"
THREADS="${THREADS:-32}"

# ----------------- LOGGING -----------------
LOG_PATH="$OUTPUT_DIR/logs"
LOG_FILE="$LOG_PATH/edna_pipeline.log"

mkdir -p "$OUTPUT_DIR" "$LOG_PATH"
touch "$LOG_FILE"

log() {
    command echo "[$(date '+%F %T')] $*" | tee -a "$LOG_FILE"
}

# Use error trap
trap 'log "ERROR: Pipeline failed at line $LINENO"; exit 1' ERR

# ------------------------------------------------------------
# Parse steps
# ------------------------------------------------------------
echo "RAW ARGS: $@" 

STEPS=()
while getopts "s:" opt; do
    case ${opt} in
        s) IFS=',' read -r -a STEPS <<< "$OPTARG" ;;
        *) echo "Usage: $0 -s <step1,step2,...>" ; exit 1 ;;
    esac
done

if [ ${#STEPS[@]} -eq 0 ]; then
    echo "Usage: $0 -s <step1,step2,...>" >&2
    exit 2
fi

log "Steps to run: ${STEPS[*]}"

# ============================================================
# STEP LOOP
# ============================================================
for STEP in "${STEPS[@]}"; do
log "Starting step: $STEP"

case "$STEP" in

# ============================================================
### QC
# ============================================================
qc)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-8}"

QC_DIR="$OUTPUT_DIR/qc"
LOG_DIR="$OUTPUT_DIR/logs/fastqc"

mkdir -p "$QC_DIR" "$LOG_DIR"

log "Starting FastQC step"
log "THREADS: $THREADS"

# -------------------- COLLECT FILES --------------------
mapfile -t R1_FILES < <(ls "$READS_DIR"/*_R1*.fastq.gz | sort)

N=${#R1_FILES[@]}
log "Total samples: $N"

# -------------------- ARRAY MODE --------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    IDX=$((SLURM_ARRAY_TASK_ID - 1))

    if [[ "$IDX" -ge "$N" ]]; then
        log "ERROR: Task ID exceeds number of samples"
        exit 1
    fi

    READS1="${R1_FILES[$IDX]}"
    READS2="${READS1/_R1/_R2}"
    SAMPLE=$(basename "$READS1" | sed -E 's/_R1.*//')

    log "Processing sample (array mode): $SAMPLE"

    fastqc \
        -t "$THREADS" \
        -o "$QC_DIR" \
        "$READS1" "$READS2" \
        &>> "$LOG_DIR/${SAMPLE}.fastqc.log"

# -------------------- NON-ARRAY FALLBACK --------------------
else

    log "No SLURM array detected → running sequentially"

    for READS1 in "${R1_FILES[@]}"; do

        READS2="${READS1/_R1/_R2}"
        SAMPLE=$(basename "$READS1" | sed -E 's/_R1.*//')

        log "Processing sample: $SAMPLE"

        fastqc \
            -t "$THREADS" \
            -o "$QC_DIR" \
            "$READS1" "$READS2" \
            &>> "$LOG_DIR/${SAMPLE}.fastqc.log"

    done

fi

log "FastQC completed"
;;

# ============================================================
### MULTIQC
# ============================================================
multiqc)

set -euo pipefail

QC_DIR="$OUTPUT_DIR/qc"
MULTIQC_DIR="$OUTPUT_DIR/multiqc"
LOG_DIR="$OUTPUT_DIR/logs/multiqc"

mkdir -p "$MULTIQC_DIR" "$LOG_DIR"

log "Starting MultiQC"

multiqc \
    "$QC_DIR" \
    -o "$MULTIQC_DIR" \
    -n "multiqc_fastqc_report.html" \
    &>> "$LOG_DIR/multiqc.log"

# -------------------- VALIDATE --------------------
if [[ ! -s "$MULTIQC_DIR/multiqc_fastqc_report.html" ]]; then
    log "ERROR: MultiQC report not generated"
    exit 1
fi

log "MultiQC report: $MULTIQC_DIR/multiqc_fastqc_report.html"
log "MultiQC completed"
;;

# ============================================================
### DEDUPLICATION (CLUMPIFY) — ARRAY-AWARE
# ============================================================
dedup)

set -euo pipefail
shopt -s nullglob

DEDUP_DIR="$OUTPUT_DIR/dedup"
LOG_DIR="$OUTPUT_DIR/logs/dedup"

mkdir -p "$DEDUP_DIR" "$LOG_DIR"

THREADS="${THREADS:-1}"

count_reads() {
    zcat "$1" | awk 'END{print NR/4}'
}

# -------------------- GET SAMPLE LIST --------------------
if [[ -n "${DEDUP_SAMPLES:-}" ]]; then

    read -ra SAMPLE_LIST <<< "$DEDUP_SAMPLES"

    R1_FILES=()

    for SAMPLE in "${SAMPLE_LIST[@]}"; do

        R1="$READS_DIR/${SAMPLE}_R1.fastq.gz"

        if [[ ! -s "$R1" ]]; then
            log "ERROR: Missing R1 file for $SAMPLE → $R1"
            exit 1
        fi

        R1_FILES+=("$R1")

    done

else

    # Default behaviour: process all samples
    mapfile -t R1_FILES < <(
        find "$READS_DIR" -maxdepth 1 \
            -type f \
            -name '*_R1*.fastq.gz' \
            | sort
    )

fi

N=${#R1_FILES[@]}

log "Number of samples selected for dedup: $N"

# -------------------- CHECK ARRAY INDEX --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: SLURM_ARRAY_TASK_ID not set → run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
READS1="${R1_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
READS2="${READS1/_R1/_R2}"
SAMPLE=$(basename "$READS1" | sed -E 's/_R1.*//')

READS1_DEDUP="$DEDUP_DIR/${SAMPLE}_dedup_R1.fastq.gz"
READS2_DEDUP="$DEDUP_DIR/${SAMPLE}_dedup_R2.fastq.gz"

log "----------------------------------------"
log "Dedup sample: $SAMPLE"
log "Task ID: ${SLURM_ARRAY_TASK_ID} / $N"

# -------------------- INPUT CHECK --------------------
[[ -s "$READS1" && -s "$READS2" ]] || {
    log "ERROR: Missing reads → exiting"
    exit 1
}

# -------------------- SKIP IF EXISTS --------------------
if [[ -s "$READS1_DEDUP" && "${FORCE:-0}" != "1" ]]; then
    log "Dedup exists → skipping"
    exit 0
fi

# -------------------- PRE-COUNTS --------------------
RAW_R1_COUNT=$(count_reads "$READS1")

# -------------------- RUN CLUMPIFY --------------------
log "Running Clumpify"

clumpify.sh \
    in1="$READS1" \
    in2="$READS2" \
    out1="$READS1_DEDUP" \
    out2="$READS2_DEDUP" \
    dedupe=t \
    optical=t \
    threads="$THREADS" \
    groups=200 \
    dupedist=12000 \
    usejni=t \
    -Xmx50g \
    &>> "$LOG_DIR/${SAMPLE}_clumpify.log"

# -------------------- VALIDATE --------------------
[[ -s "$READS1_DEDUP" && -s "$READS2_DEDUP" ]] || {
    log "ERROR: Dedup failed"
    exit 1
}

# -------------------- POST-COUNTS --------------------
DEDUP_R1_COUNT=$(count_reads "$READS1_DEDUP")

log "After dedup: R1=${DEDUP_R1_COUNT}"

if [[ "$RAW_R1_COUNT" -gt 0 ]]; then
    PCT_REMOVED=$(awk -v raw="$RAW_R1_COUNT" -v dedup="$DEDUP_R1_COUNT" \
        'BEGIN {printf "%.2f", (1 - dedup/raw)*100}')
    log "Duplicate removal: ${PCT_REMOVED}%"
fi

log "Finished sample: $SAMPLE"
;;

# ============================================================
### TRIMMING (FASTP — ARRAY-AWARE)
# ============================================================
trim)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-1}"

module --force purge
module load cluster/genius/batch
module load fastp/1.0.1-GCC-13.3.0

DEDUP_DIR="$OUTPUT_DIR/dedup"
TRIM_DIR="$OUTPUT_DIR/trimmed"
LOG_DIR="$OUTPUT_DIR/logs/fastp"

mkdir -p "$TRIM_DIR" "$LOG_DIR"

log "Starting trimming"
log "Input:  $DEDUP_DIR"
log "Output: $TRIM_DIR"
log "Threads: $THREADS"

count_reads() {
    zcat "$1" | awk 'END{print NR/4}'
}

# -------------------- GET SAMPLE LIST --------------------
if [[ -n "${TRIM_SAMPLES:-}" ]]; then

    # Use only explicitly requested samples
    read -ra SAMPLE_LIST <<< "$TRIM_SAMPLES"

    R1_FILES=()

    for SAMPLE in "${SAMPLE_LIST[@]}"; do

        R1="$DEDUP_DIR/${SAMPLE}_dedup_R1.fastq.gz"
        R2="$DEDUP_DIR/${SAMPLE}_dedup_R2.fastq.gz"

        # Check both paired files before adding sample
        if [[ ! -s "$R1" ]]; then
            log "ERROR: Missing R1 for $SAMPLE → $R1"
            exit 1
        fi

        if [[ ! -s "$R2" ]]; then
            log "ERROR: Missing R2 for $SAMPLE → $R2"
            exit 1
        fi

        R1_FILES+=("$R1")
    done

else

    # Default behaviour: use all deduplicated samples
    mapfile -t R1_FILES < <(
        find "$DEDUP_DIR" \
            -maxdepth 1 \
            -type f \
            -name '*_dedup_R1.fastq.gz' \
            | sort
    )

fi

N=${#R1_FILES[@]}

log "Number of samples selected for trimming: $N"

# -------------------- CHECK ARRAY --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: SLURM_ARRAY_TASK_ID not set → use array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
READS1="${R1_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
READS2="${READS1/_R1.fastq.gz/_R2.fastq.gz}"
SAMPLE=$(basename "$READS1" | sed -E 's/_dedup_R1.*//')

READS1_TRIM="$TRIM_DIR/${SAMPLE}_trimmed_R1.fastq.gz"
READS2_TRIM="$TRIM_DIR/${SAMPLE}_trimmed_R2.fastq.gz"

HTML="$LOG_DIR/${SAMPLE}.html"
JSON="$LOG_DIR/${SAMPLE}.json"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing $SAMPLE (Task ${SLURM_ARRAY_TASK_ID}/$N)"

# -------------------- INPUT CHECK --------------------
[[ -s "$READS1" && -s "$READS2" ]] || {
    log "ERROR: Missing input → exiting"
    exit 1
}

# -------------------- SKIP --------------------
if [[ -s "$READS1_TRIM" && "${FORCE:-0}" != "1" ]]; then
    log "Trim exists → skipping"
    exit 0
fi

# -------------------- PRE-COUNT --------------------
BEFORE=$(count_reads "$READS1")

# -------------------- RUN FASTP --------------------
log "Running fastp"

fastp \
    -i "$READS1" \
    -I "$READS2" \
    -o "$READS1_TRIM" \
    -O "$READS2_TRIM" \
    --detect_adapter_for_pe \
    --thread "$THREADS" \
    --overrepresentation_analysis \
    --trim_poly_g \
    --poly_g_min_len 30 \
    --correction \
    --cut_front \
    --cut_tail \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html "$HTML" \
    --json "$JSON" \
    &>> "$LOGFILE"

# -------------------- VALIDATE --------------------
[[ -s "$READS1_TRIM" && -s "$READS2_TRIM" ]] || {
    log "ERROR: fastp failed"
    exit 1
}

# -------------------- POST-COUNT --------------------
AFTER=$(count_reads "$READS1_TRIM")

log "Reads: ${BEFORE} → ${AFTER}"

if [[ "$BEFORE" -gt 0 ]]; then
    PCT=$(awk -v b="$BEFORE" -v a="$AFTER" \
        'BEGIN {printf "%.2f", (1 - a/b)*100}')
    log "Trimmed: ${PCT}%"
fi

log "Finished sample: $SAMPLE"
;;

# ============================================================
### MERGE READS (INTERLEAVED — ARRAY-AWARE)
# ============================================================
merge_reads)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-1}"

TRIM_DIR="$OUTPUT_DIR/trimmed"
MERGE_DIR="$OUTPUT_DIR/merged_reads"
LOG_DIR="$OUTPUT_DIR/logs/merge_reads"

mkdir -p "$MERGE_DIR" "$LOG_DIR"

log "Starting read merging (interleaved)"

# -------------------- TOOL CHECK --------------------
command -v reformat.sh >/dev/null || {
    log "ERROR: BBTools (reformat.sh) not found"
    exit 1
}

# -------------------- GET SAMPLE LIST --------------------
mapfile -t R1_FILES < <(ls "$TRIM_DIR"/*_trimmed_R1.fastq.gz | sort)
N=${#R1_FILES[@]}

[[ "$N" -gt 0 ]] || {
    log "ERROR: No trimmed reads found"
    exit 1
}

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: SLURM_ARRAY_TASK_ID not set → use array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
READS1="${R1_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
READS2="${READS1/_R1.fastq.gz/_R2.fastq.gz}"
SAMPLE=$(basename "$READS1" | sed -E 's/_trimmed_R1.*//')

MERGED="$MERGE_DIR/${SAMPLE}.interleaved.fastq.gz"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing: $SAMPLE"
log "Input: $READS1 $READS2"
log "Output: $MERGED"

# -------------------- INPUT CHECK --------------------
[[ -s "$READS1" && -s "$READS2" ]] || {
    log "ERROR: Missing trimmed reads"
    exit 1
}

# -------------------- SKIP --------------------
if [[ -s "$MERGED" && "${FORCE:-0}" != "1" ]]; then
    log "Merged reads exist → skipping"
    exit 0
fi

# -------------------- RUN MERGE --------------------
log "Merging reads (interleaved)"

reformat.sh \
    in1="$READS1" \
    in2="$READS2" \
    out="$MERGED" \
    threads="$THREADS" \
    &>> "$LOGFILE"

# -------------------- VALIDATE --------------------
[[ -s "$MERGED" ]] || {
    log "ERROR: Merging failed"
    exit 1
}

log "Finished merging: $SAMPLE"
;;

# ============================================================
### POST-TRIM QC (FASTQC — ARRAY-AWARE)
# ============================================================
postqc)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-4}"

TRIM_DIR="$OUTPUT_DIR/trimmed"
POSTQC_DIR="$OUTPUT_DIR/postqc"
LOG_DIR="$OUTPUT_DIR/logs/postqc"

mkdir -p "$POSTQC_DIR" "$LOG_DIR"

log "Starting post-trim FastQC"
log "Input:  $TRIM_DIR"
log "Output: $POSTQC_DIR"
log "Threads: $THREADS"

# -------------------- MODULES --------------------
module --force purge
module load cluster/genius/batch
module load FastQC

# -------------------- GET SAMPLE LIST --------------------
mapfile -t R1_FILES < <(ls "$TRIM_DIR"/*_trimmed_R1.fastq.gz | sort)
N=${#R1_FILES[@]}

log "Total samples: $N"

# -------------------- ARRAY MODE --------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    IDX=$((SLURM_ARRAY_TASK_ID - 1))

    if [[ "$IDX" -ge "$N" ]]; then
        log "ERROR: Task ID exceeds number of samples"
        exit 1
    fi

    READS1="${R1_FILES[$IDX]}"
    READS2="${READS1/_R1.fastq.gz/_R2.fastq.gz}"
    SAMPLE=$(basename "$READS1" | sed -E 's/_trimmed_R1.*//')

    log "Processing sample (array mode): $SAMPLE"

    fastqc \
        -t "$THREADS" \
        -o "$POSTQC_DIR" \
        "$READS1" "$READS2" \
        &>> "$LOG_DIR/${SAMPLE}.postqc.log"

# -------------------- NON-ARRAY FALLBACK --------------------
else

    log "No SLURM array detected → running sequentially"

    for READS1 in "${R1_FILES[@]}"; do

        READS2="${READS1/_R1.fastq.gz/_R2.fastq.gz}"
        SAMPLE=$(basename "$READS1" | sed -E 's/_trimmed_R1.*//')

        log "Processing sample: $SAMPLE"

        fastqc \
            -t "$THREADS" \
            -o "$POSTQC_DIR" \
            "$READS1" "$READS2" \
            &>> "$LOG_DIR/${SAMPLE}.postqc.log"

    done

fi

log "Post-trim FastQC completed"
;;

# ============================================================
### KRAKEN TAXONOMY — STANDARD DB (ARRAY-AWARE)
# ============================================================
kraken_standard)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-1}"

# -------------------- FIXED: removed space after :- --------------------
KRAKEN_DB="${KRAKEN_DB:-/staging/leuven/stg_00100/Rodgee/databases/kraken2_standard_20260226/db}"

KRAKEN_OUT="$OUTPUT_DIR/kraken2_standard"
LOG_DIR="$OUTPUT_DIR/logs/kraken2_standard"

mkdir -p "$KRAKEN_OUT" "$LOG_DIR"

log "Starting Kraken STANDARD step"
log "DB: $KRAKEN_DB"
log "THREADS: $THREADS"

# -------------------- CHECK DB --------------------
[[ -d "$KRAKEN_DB" ]] || {
    log "ERROR: Kraken STANDARD DB not found: $KRAKEN_DB"
    exit 1
}

# -------------------- COLLECT FILES --------------------
R1_FILES=("$OUTPUT_DIR"/trimmed/*_R1.fastq.gz)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    log "ERROR: No input FASTQ files found in $OUTPUT_DIR/trimmed"
    exit 1
fi

# -------------------- SORT SAFELY --------------------
IFS=$'\n' read -r -d '' -a R1_FILES < <(printf "%s\n" "${R1_FILES[@]}" | sort && printf '\0')

N=${#R1_FILES[@]}
log "Total samples: $N"

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
R1="${R1_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$R1" _R1.fastq.gz)
R2="$OUTPUT_DIR/trimmed/${SAMPLE}_R2.fastq.gz"

REPORT="$KRAKEN_OUT/${SAMPLE}.report"
OUT="$KRAKEN_OUT/${SAMPLE}.out"

CLASS1="$KRAKEN_OUT/${SAMPLE}_classified_1.fastq"
CLASS2="$KRAKEN_OUT/${SAMPLE}_classified_2.fastq"
UNCLASS1="$KRAKEN_OUT/${SAMPLE}_unclassified_1.fastq"
UNCLASS2="$KRAKEN_OUT/${SAMPLE}_unclassified_2.fastq"

log "----------------------------------------"
log "Processing $SAMPLE (STANDARD DB)"

# -------------------- INPUT CHECK --------------------
[[ -s "$R1" && -s "$R2" ]] || {
    log "ERROR: Missing input reads for $SAMPLE"
    exit 1
}

# -------------------- SKIP --------------------
if [[ -s "$REPORT" && "${FORCE:-0}" != "1" ]]; then
    log "Kraken STANDARD already done → skipping"
    exit 0
fi

# -------------------- RUN --------------------
log "Running Kraken2 (STANDARD)"

kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired \
    --confidence 0.1 \
    --use-names \
    --report "$REPORT" \
    --output "$OUT" \
    --classified-out "$KRAKEN_OUT/${SAMPLE}_classified#.fastq" \
    --unclassified-out "$KRAKEN_OUT/${SAMPLE}_unclassified#.fastq" \
    "$R1" "$R2" \
    2>&1 | tee "$LOG_DIR/${SAMPLE}.kraken.log"

# -------------------- VALIDATE --------------------
[[ -s "$REPORT" ]] || {
    log "ERROR: Kraken STANDARD failed for $SAMPLE"
    exit 1
}

# -------------------- COMPRESS --------------------
log "Compressing outputs"

for f in "$CLASS1" "$CLASS2" "$UNCLASS1" "$UNCLASS2"; do
    [[ -s "$f" ]] && gzip -f "$f"
done

log "Finished STANDARD: $SAMPLE"
;;

# ============================================================
# PER-SAMPLE ASSEMBLY (ARRAY OR SINGLE-SAMPLE)
# ============================================================
assembly_sample)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-8}"
MINLEN="${MINLEN:-1000}"

TRIM_DIR="${TRIM_DIR:-$OUTPUT_DIR/trimmed}"
ASM_DIR="$OUTPUT_DIR/assembly"
LOG_DIR="$OUTPUT_DIR/logs/assembly"

mkdir -p "$ASM_DIR" "$LOG_DIR"

log "Starting assembly step"
log "THREADS: $THREADS"
log "MINLEN: $MINLEN"
log "TRIM_DIR: $TRIM_DIR"

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    R1="$TRIM_DIR/${SAMPLE}_trimmed_R1.fastq.gz"
    R2="$TRIM_DIR/${SAMPLE}_trimmed_R2.fastq.gz"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    mapfile -t R1_FILES < <(
        find "$TRIM_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*_trimmed_R1.fastq.gz" \
            | sort
    )

    N=${#R1_FILES[@]}

    log "Total samples detected: $N"

    # -------------------- CHECK FILE LIST --------------------
    if [[ "$N" -eq 0 ]]; then
        log "ERROR: No trimmed R1 files found in:"
        log "$TRIM_DIR"
        exit 1
    fi

    # -------------------- CHECK ARRAY --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    R1="${R1_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

    SAMPLE=$(basename "$R1" | sed -E 's/_trimmed_R1.*//')

fi

# ============================================================
# INPUT CHECK
# ============================================================

[[ -s "$R1" ]] || {
    log "ERROR: Missing or empty R1:"
    log "$R1"
    exit 1
}

[[ -s "$R2" ]] || {
    log "ERROR: Missing or empty R2:"
    log "$R2"
    exit 1
}

# ============================================================
# OUTPUT FILES
# ============================================================

OUTDIR="$ASM_DIR/$SAMPLE"

CONTIGS="$OUTDIR/final.contigs.fa"
FINAL="$OUTDIR/${SAMPLE}.contigs.fa"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "R1: $R1"
log "R2: $R2"
log "Output directory: $OUTDIR"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# SKIP IF DONE
# ============================================================

if [[ -s "$FINAL" && "${FORCE:-0}" != "1" ]]; then
    log "Assembly already exists -> skipping"
    log "$FINAL"
    exit 0
fi

# ============================================================
# CLEAN OLD/INCOMPLETE RUN
# ============================================================

if [[ -d "$OUTDIR" ]]; then
    log "Removing existing assembly directory:"
    log "$OUTDIR"
    rm -rf "$OUTDIR"
fi

# Do not mkdir OUTDIR.
# MEGAHIT requires the output directory not to already exist.

# ============================================================
# RUN MEGAHIT
# ============================================================

log "Running MEGAHIT"

megahit \
    -1 "$R1" \
    -2 "$R2" \
    -t "$THREADS" \
    --min-contig-len "$MINLEN" \
    --mem-flag 0 \
    -o "$OUTDIR" \
    &>> "$LOG_DIR/${SAMPLE}.log"

# ============================================================
# VALIDATE ASSEMBLY
# ============================================================

[[ -s "$CONTIGS" ]] || {
    log "ERROR: Assembly failed or final.contigs.fa is empty"
    log "Check:"
    log "$LOG_DIR/${SAMPLE}.log"
    exit 1
}

# ============================================================
# RENAME CONTIG HEADERS
# ============================================================

TMP="$OUTDIR/${SAMPLE}.renamed.fa"

awk -v s="$SAMPLE" '
/^>/ {
    sub(/^>/, "", $0)
    print ">" s "_" $0
    next
}
{
    print
}
' "$CONTIGS" > "$TMP"

# Rename the prefixed-contig file to the final sample-specific name
mv "$TMP" "$FINAL"

# Replace MEGAHIT final.contigs.fa with a symlink to the
# sample-specific renamed assembly
ln -sfn "$(basename "$FINAL")" "$CONTIGS"

# ============================================================
# FINAL CHECK
# ============================================================

[[ -s "$FINAL" ]] || {
    log "ERROR: Final renamed assembly not created"
    exit 1
}

log "Assembly completed"
log "Final assembly: $FINAL"

;;

# ============================================================
# PER-SAMPLE MAPPING (ARRAY OR SINGLE-SAMPLE)
# ============================================================
mapping)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-8}"

MAP_DIR="$OUTPUT_DIR/mapping"
LOG_DIR="$OUTPUT_DIR/logs/mapping"

mkdir -p "$MAP_DIR" "$LOG_DIR"

log "Starting mapping step (BBMap)"
log "THREADS: $THREADS"

# ============================================================
# TMPDIR VALIDATION
# ============================================================

if [[ -z "${TMPDIR:-}" || ! -d "$TMPDIR" ]]; then
    log "ERROR: TMPDIR not set or not a directory"
    exit 1
fi

# Test whether TMPDIR is writable
touch "$TMPDIR/.tmp_test_${SLURM_JOB_ID}" 2>/dev/null || {
    log "ERROR: TMPDIR not writable: $TMPDIR"
    exit 1
}

rm -f "$TMPDIR/.tmp_test_${SLURM_JOB_ID}"

# Check available temporary space
FREE_GB=$(df -BG "$TMPDIR" | awk 'NR==2{gsub("G","",$4); print $4}')
MIN_GB="${TMP_MIN_GB:-20}"

if (( FREE_GB < MIN_GB )); then
    log "ERROR: TMPDIR has only ${FREE_GB}G free (< ${MIN_GB}G required)"
    exit 1
fi

# Per-job working directory
WORKDIR="$TMPDIR/map_${SLURM_JOB_ID}"

mkdir -p "$WORKDIR"

log "Using TMPDIR: $TMPDIR"
log "WORKDIR: $WORKDIR"
log "TMP free space: ${FREE_GB}G"

# ============================================================
# TOOL CHECK
# ============================================================

command -v bbmap.sh >/dev/null || {
    log "ERROR: bbmap.sh not found"
    exit 1
}

command -v samtools >/dev/null || {
    log "ERROR: samtools not found"
    exit 1
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    ASM="$OUTPUT_DIR/assembly/$SAMPLE/${SAMPLE}.contigs.fa"

    [[ -s "$ASM" ]] || {
        log "ERROR: Missing assembly:"
        log "$ASM"
        exit 1
    }

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$OUTPUT_DIR/assembly" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.fa" \
            | sort
    )

    N=${#ASM_FILES[@]}

    log "Total assemblies detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No assemblies found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$ASM" .contigs.fa)

fi

# ============================================================
# SAMPLE FILES
# ============================================================

R1="$OUTPUT_DIR/trimmed/${SAMPLE}_trimmed_R1.fastq.gz"
R2="$OUTPUT_DIR/trimmed/${SAMPLE}_trimmed_R2.fastq.gz"

BAM="$MAP_DIR/${SAMPLE}.bam"
BAI="$BAM.bai"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "Assembly: $ASM"
log "R1: $R1"
log "R2: $R2"
log "BAM: $BAM"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# INPUT CHECK
# ============================================================

[[ -s "$ASM" ]] || {
    log "ERROR: Missing assembly:"
    log "$ASM"
    exit 1
}

[[ -s "$R1" ]] || {
    log "ERROR: Missing R1:"
    log "$R1"
    exit 1
}

[[ -s "$R2" ]] || {
    log "ERROR: Missing R2:"
    log "$R2"
    exit 1
}

# ============================================================
# SKIP IF DONE
# ============================================================

if [[ -s "$BAM" && -s "$BAI" && "${FORCE:-0}" != "1" ]]; then
    log "Mapping already completed -> skipping"
    log "$BAM"
    exit 0
fi

# ============================================================
# CLEAN PARTIAL OUTPUT
# ============================================================

if [[ "${FORCE:-0}" == "1" ]]; then
    log "FORCE=1 -> removing previous BAM and index"
    rm -f "$BAM" "$BAI"
fi

# ============================================================
# TMP FILES
# ============================================================

TMP_PREFIX="$WORKDIR/$SAMPLE"

TMP_REF="${TMP_PREFIX}.fa"
TMP_R1="${TMP_PREFIX}_R1.fastq.gz"
TMP_R2="${TMP_PREFIX}_R2.fastq.gz"

TMP_BAM="${TMP_PREFIX}.bam"

log "Copying inputs to WORKDIR"

cp "$ASM" "$TMP_REF"
cp "$R1" "$TMP_R1"
cp "$R2" "$TMP_R2"

# ============================================================
# THREAD SPLIT
# ============================================================

BBMAP_THREADS=$((THREADS - 4))
SORT_THREADS=4

(( BBMAP_THREADS < 1 )) && BBMAP_THREADS=1

log "BBMap threads: $BBMAP_THREADS"
log "samtools sort threads: $SORT_THREADS"

# ============================================================
# JAVA MEMORY
# ============================================================

MEM_MB="${SLURM_MEM_PER_NODE:-0}"

if (( MEM_MB > 0 )); then

    MEM_GB=$((MEM_MB / 1024))

    # Give BBMap 80% of allocated memory
    JAVA_MEM=$((MEM_GB * 80 / 100))

else

    JAVA_MEM=16

fi

log "Java memory: ${JAVA_MEM}g"

# ============================================================
# RUN BBMAP
# ============================================================

log "Running BBMap"

bbmap.sh \
    ref="$TMP_REF" \
    in1="$TMP_R1" \
    in2="$TMP_R2" \
    threads="$BBMAP_THREADS" \
    minid=0.95 \
    path="$WORKDIR/bbmap_${SAMPLE}" \
    -Xmx${JAVA_MEM}g \
    out=stdout.sam \
    2>> "$LOG_DIR/${SAMPLE}.log" \
| samtools view \
    -b \
    - \
| samtools sort \
    -@ "$SORT_THREADS" \
    -T "$WORKDIR/${SAMPLE}" \
    -o "$TMP_BAM" \
    2>> "$LOG_DIR/${SAMPLE}.log"

# ============================================================
# VALIDATE TEMP BAM
# ============================================================

[[ -s "$TMP_BAM" ]] || {
    log "ERROR: BAM not created"
    exit 1
}

# Quick BAM integrity check
samtools quickcheck "$TMP_BAM" || {
    log "ERROR: BAM failed samtools quickcheck"
    exit 1
}

# ============================================================
# COPY BACK
# ============================================================

log "Copying BAM to output directory"

cp "$TMP_BAM" "$BAM"

# ============================================================
# INDEX BAM
# ============================================================

log "Indexing BAM"

samtools index \
    -@ "$SORT_THREADS" \
    "$BAM"

[[ -s "$BAI" ]] || {
    log "ERROR: BAM index not created"
    exit 1
}

# ============================================================
# FINAL VALIDATION
# ============================================================

samtools quickcheck "$BAM" || {
    log "ERROR: Final BAM failed samtools quickcheck"
    exit 1
}

# ============================================================
# CLEANUP
# ============================================================

if [[ "${DEBUG:-0}" != "1" ]]; then

    log "Removing temporary working directory"

    rm -rf "$WORKDIR"

else

    log "DEBUG=1 -> keeping WORKDIR:"
    log "$WORKDIR"

fi

log "Finished mapping for $SAMPLE"

;;

# ============================================================
# QUAST (ASSEMBLY QC — ARRAY OR SINGLE-SAMPLE)
# ============================================================
quast)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-8}"

ASM_DIR="$OUTPUT_DIR/assembly"
QUAST_DIR="$OUTPUT_DIR/quast"
LOG_DIR="$OUTPUT_DIR/logs/quast"

mkdir -p "$QUAST_DIR" "$LOG_DIR"

log "Starting QUAST step"
log "THREADS: $THREADS"
log "ASM_DIR: $ASM_DIR"
log "QUAST_DIR: $QUAST_DIR"

# ============================================================
# TOOL CHECK
# ============================================================

command -v quast.py >/dev/null 2>&1 || {
    log "ERROR: quast.py not found"
    log "Make sure QUAST is loaded/activated in the SLURM job"
    exit 1
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    ASM="$ASM_DIR/$SAMPLE/${SAMPLE}.contigs.fa"

    [[ -s "$ASM" ]] || {
        log "ERROR: Assembly not found:"
        log "$ASM"
        exit 1
    }

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$ASM_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.fa" \
            | sort
    )

    N=${#ASM_FILES[@]}

    log "Total valid assemblies detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No valid assemblies found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$ASM" .contigs.fa)

fi

# ============================================================
# OUTPUT
# ============================================================

OUTDIR="$QUAST_DIR/$SAMPLE"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "Assembly: $ASM"
log "Output: $OUTDIR"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# INPUT CHECK
# ============================================================

[[ -s "$ASM" ]] || {
    log "ERROR: Missing or empty assembly:"
    log "$ASM"
    exit 1
}

# ============================================================
# SKIP IF DONE
# ============================================================

if [[ -s "$OUTDIR/report.txt" && "${FORCE:-0}" != "1" ]]; then
    log "QUAST already completed -> skipping"
    log "$OUTDIR/report.txt"
    exit 0
fi

# ============================================================
# CLEAN OLD / PARTIAL RUN
# ============================================================

if [[ -d "$OUTDIR" ]]; then
    log "Removing previous QUAST directory:"
    log "$OUTDIR"
    rm -rf "$OUTDIR"
fi

# ============================================================
# RUN QUAST
# ============================================================

log "Running QUAST"

quast.py \
    "$ASM" \
    -o "$OUTDIR" \
    -t "$THREADS" \
    --min-contig 1000 \
    &>> "$LOG_DIR/${SAMPLE}.log"

# ============================================================
# VALIDATE
# ============================================================

if [[ ! -s "$OUTDIR/report.txt" ]]; then
    log "ERROR: QUAST failed for $SAMPLE"
    log "Check log:"
    log "$LOG_DIR/${SAMPLE}.log"
    exit 1
fi

# ============================================================
# SUMMARY
# ============================================================

log "QUAST completed for $SAMPLE"
log "Report: $OUTDIR/report.txt"

;;

# ============================================================
# QUAST SUMMARY (MERGE ALL SAMPLES INTO TSV)
# ============================================================
quast_summary)

set -euo pipefail
shopt -s nullglob

QUAST_DIR="$OUTPUT_DIR/quast"
META="/staging/leuven/stg_00100/Rodgee/caves_data/metadata/biological_metadata.tsv"
OUT="$QUAST_DIR/quast_summary.tsv"

mkdir -p "$QUAST_DIR"

log "Creating QUAST summary table"

# ============================================================
# INPUT CHECKS
# ============================================================

[[ -d "$QUAST_DIR" ]] || {
    log "ERROR: QUAST directory not found:"
    log "$QUAST_DIR"
    exit 1
}

[[ -s "$META" ]] || {
    log "ERROR: Metadata file not found or empty:"
    log "$META"
    exit 1
}

# ============================================================
# HEADER
# ============================================================

printf "sample\tsample_code\tincubation\tn_contigs\ttotal_length\tN50\tGC\n" > "$OUT"

# ============================================================
# FIND REPORTS
# ============================================================

mapfile -t REPORTS < <(
    find "$QUAST_DIR" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name "report.txt" \
        | sort
)

NREPORTS=${#REPORTS[@]}

log "QUAST reports detected: $NREPORTS"

if (( NREPORTS == 0 )); then
    log "ERROR: No QUAST report.txt files found"
    exit 1
fi

# ============================================================
# LOOP THROUGH SAMPLES
# ============================================================

NWRITTEN=0

for REPORT in "${REPORTS[@]}"; do

    SAMPLE=$(basename "$(dirname "$REPORT")")

    log "Processing: $SAMPLE"

    # --------------------------------------------------------
    # METADATA
    # --------------------------------------------------------

    INC="NA"
    CODE="NA"

    META_LINE=$(awk -F'\t' -v s="$SAMPLE" '
        $1 == s {
            print $2 "\t" $3
            exit
        }
    ' "$META")

    if [[ -n "$META_LINE" ]]; then
        IFS=$'\t' read -r INC CODE <<< "$META_LINE"
    else
        log "WARNING: No metadata found for $SAMPLE"
    fi

    # --------------------------------------------------------
    # QUAST METRICS
    # --------------------------------------------------------

    ncontigs=$(awk '
        /^# contigs \(>= 0 bp\)/ {
            print $NF
            exit
        }
    ' "$REPORT")

    totlen=$(awk '
        /^Total length \(>= 0 bp\)/ {
            print $NF
            exit
        }
    ' "$REPORT")

    n50=$(awk '
        /^N50[[:space:]]/ {
            print $NF
            exit
        }
    ' "$REPORT")

    gc=$(awk '
        /^GC \(%\)/ {
            print $NF
            exit
        }
    ' "$REPORT")

    # --------------------------------------------------------
    # DEFAULT TO NA IF METRIC IS ABSENT
    # --------------------------------------------------------

    ncontigs="${ncontigs:-NA}"
    totlen="${totlen:-NA}"
    n50="${n50:-NA}"
    gc="${gc:-NA}"

    # --------------------------------------------------------
    # WRITE OUTPUT
    # --------------------------------------------------------

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE" \
        "$CODE" \
        "$INC" \
        "$ncontigs" \
        "$totlen" \
        "$n50" \
        "$gc" \
        >> "$OUT"

    ((NWRITTEN+=1))

done

# ============================================================
# SUMMARY
# ============================================================

log "Created QUAST summary:"
log "$OUT"

log "Samples written: $NWRITTEN"

;;

# ============================================================
# PER-SAMPLE COVERAGE (COVERM — ARRAY OR SINGLE-SAMPLE)
# ============================================================
coverage_coverm)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-8}"

ASM_DIR="$OUTPUT_DIR/assembly"
MAP_DIR="$OUTPUT_DIR/mapping"
COV_DIR="$OUTPUT_DIR/coverage_coverm"
LOG_DIR="$OUTPUT_DIR/logs/coverage_coverm"

mkdir -p "$COV_DIR" "$LOG_DIR"

log "Starting CoverM coverage step"
log "THREADS: $THREADS"
log "ASM_DIR: $ASM_DIR"
log "MAP_DIR: $MAP_DIR"
log "COV_DIR: $COV_DIR"

# ============================================================
# TOOL CHECK
# ============================================================

command -v coverm >/dev/null 2>&1 || {
    log "ERROR: coverm not found"
    log "Make sure the CoverM environment is activated in SLURM"
    exit 1
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    ASM="$ASM_DIR/$SAMPLE/${SAMPLE}.contigs.fa"

    [[ -s "$ASM" ]] || {
        log "ERROR: Assembly not found:"
        log "$ASM"
        exit 1
    }

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$ASM_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.fa" \
            | sort
    )

    N=${#ASM_FILES[@]}

    log "Total assemblies detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No valid assemblies found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$ASM" .contigs.fa)

fi

# ============================================================
# SAMPLE FILES
# ============================================================

BAM="$MAP_DIR/${SAMPLE}.bam"

RAW="$COV_DIR/${SAMPLE}.coverage.raw.tsv"
OUT="$COV_DIR/${SAMPLE}.coverage.tsv"

LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "Assembly: $ASM"
log "BAM: $BAM"
log "Raw coverage: $RAW"
log "Final coverage: $OUT"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# INPUT CHECK
# ============================================================

[[ -s "$ASM" ]] || {
    log "ERROR: Missing assembly:"
    log "$ASM"
    exit 1
}

[[ -s "$BAM" ]] || {
    log "ERROR: Missing BAM:"
    log "$BAM"
    exit 1
}

# Check BAM integrity before using it
samtools quickcheck "$BAM" 2>/dev/null || {
    log "ERROR: BAM failed integrity check:"
    log "$BAM"
    exit 1
}

# ============================================================
# STEP 1: RUN COVERM
# ============================================================

if [[ -s "$RAW" && "${FORCE:-0}" != "1" ]]; then

    log "Raw coverage already exists -> skipping CoverM"

else

    # Remove partial/old raw result
    rm -f "$RAW"

    log "Running CoverM"

    coverm contig \
        --bam-files "$BAM" \
        --threads "$THREADS" \
        --methods mean covered_fraction length \
        --min-read-percent-identity 0.95 \
        --min-read-aligned-percent 75 \
        --output-format dense \
        > "$RAW" \
        2> "$LOGFILE"

    # -------------------- VALIDATE RAW --------------------
    [[ -s "$RAW" ]] || {
        log "ERROR: CoverM failed or produced empty output"
        log "Check log:"
        log "$LOGFILE"
        exit 1
    }

fi

# ============================================================
# STEP 2: COMPUTE WEIGHTED COVERAGE
# ============================================================

if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then

    log "Final coverage already exists -> skipping weighting"
    exit 0

fi

# Remove partial previous output
rm -f "$OUT"

log "Computing coverage_weighted"

# Save CoverM header to log for debugging
log "CoverM header:"
head -n 1 "$RAW" | tee -a "$LOGFILE"

# ============================================================
# ROBUST COVERM OUTPUT PARSER
# ============================================================

awk '
BEGIN {
    OFS="\t"
}

NR==1 {
    print "contig", "mean_coverage", "covered_fraction", "length", "coverage_weighted"
    next
}

{
    contig=$1

    # Reset numeric-field counter for every row
    n=0

    # Collect numeric fields from CoverM output
    for (i=1; i<=NF; i++) {

        if ($i ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
            nums[++n]=$i
        }
    }

    # Expected final numeric fields:
    # mean, covered_fraction, length
    if (n < 3)
        next

    cov=nums[n-2]
    frac=nums[n-1]
    len=nums[n]

    weighted=cov*len

    print contig, cov, frac, len, weighted
}
' "$RAW" > "$OUT"

# ============================================================
# VALIDATE FINAL OUTPUT
# ============================================================

[[ -s "$OUT" ]] || {
    log "ERROR: Final coverage output not created"
    exit 1
}

NLINES=$(wc -l < "$OUT")

if (( NLINES <= 1 )); then
    log "WARNING: Coverage output contains header only"
fi

# ============================================================
# SUMMARY
# ============================================================

log "Coverage preview:"
head -n 5 "$OUT" | tee -a "$LOGFILE"

log "Number of output lines: $NLINES"
log "Finished coverage for $SAMPLE"

;;

# ============================================================
# COVERAGE SUMMARY PER SAMPLE
# ============================================================
coverage_summary_sample)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

COV_DIR="$OUTPUT_DIR/coverage_coverm"
OUT="$OUTPUT_DIR/coverage_coverm/coverage_per_sample.tsv"

log "Summarizing coverage per sample (weighted)"

echo -e "sample\ttotal_cov_weighted\ttotal_length\tmean_coverage" > "$OUT"

for f in "$COV_DIR"/*.coverage.tsv; do

    sample=$(basename "$f" .coverage.tsv)

    awk '
    NR>1 {
        cov += $5      # coverage_weighted
        len += $4      # length
    }
    END {
        if (len > 0)
            mean = cov / len
        else
            mean = 0

        printf "%.6f\t%d\t%.6f\n", cov, len, mean
    }
    ' "$f" | awk -v s="$sample" '{print s "\t" $0}' >> "$OUT"

done

log "Created: $OUT"

;;

# ============================================================
# PER-SAMPLE COVERAGE (SAMTOOLS -aa — ARRAY OR SINGLE-SAMPLE)
# ============================================================
coverage_samtools)

set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-8}"

MAP_DIR="$OUTPUT_DIR/mapping"
COV_DIR="$OUTPUT_DIR/coverage_samtools"
LOG_DIR="$OUTPUT_DIR/logs/coverage_samtools"

mkdir -p "$COV_DIR" "$LOG_DIR"

log "Starting samtools coverage step"
log "THREADS: $THREADS"
log "MAP_DIR: $MAP_DIR"
log "COV_DIR: $COV_DIR"

# ============================================================
# TOOL CHECK
# ============================================================

command -v samtools >/dev/null 2>&1 || {
    log "ERROR: samtools not found"
    exit 1
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    BAM="$MAP_DIR/${SAMPLE}.bam"

    [[ -s "$BAM" ]] || {
        log "ERROR: BAM not found:"
        log "$BAM"
        exit 1
    }

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    mapfile -t BAM_FILES < <(
        find "$MAP_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.bam" \
            -size +0c \
            | sort
    )

    N=${#BAM_FILES[@]}

    log "Total BAM files detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No BAM files found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    BAM="${BAM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$BAM" .bam)

fi

# ============================================================
# SAMPLE FILES
# ============================================================

OUT="$COV_DIR/${SAMPLE}.samtools.tsv"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "BAM: $BAM"
log "Output: $OUT"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# INPUT CHECK
# ============================================================

[[ -s "$BAM" ]] || {
    log "ERROR: Missing or empty BAM:"
    log "$BAM"
    exit 1
}

# Check BAM integrity
samtools quickcheck "$BAM" || {
    log "ERROR: BAM failed samtools quickcheck:"
    log "$BAM"
    exit 1
}

# ============================================================
# SKIP IF DONE
# ============================================================

if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
    log "Samtools coverage already completed -> skipping"
    log "$OUT"
    exit 0
fi

# ============================================================
# CLEAN PARTIAL OUTPUT
# ============================================================

rm -f "$OUT"

# ============================================================
# RUN SAMTOOLS DEPTH
# ============================================================

log "Running samtools depth (-aa)"

samtools depth \
    -aa \
    "$BAM" \
    2>> "$LOGFILE" \
| awk -v sample="$SAMPLE" '
BEGIN {
    sum = 0
    total = 0
    covered = 0
}

{
    depth = $3

    total++
    sum += depth

    if (depth > 0) {
        covered++
    }
}

END {

    if (total > 0) {
        mean = sum / total
        breadth = covered / total
    }
    else {
        mean = 0
        breadth = 0
    }

    print "sample\tmean_coverage\tbreadth"
    printf "%s\t%.6f\t%.6f\n", sample, mean, breadth
}
' > "$OUT"

# ============================================================
# VALIDATE
# ============================================================

[[ -s "$OUT" ]] || {
    log "ERROR: samtools coverage output not created"
    exit 1
}

NLINES=$(wc -l < "$OUT")

if (( NLINES < 2 )); then
    log "ERROR: Coverage output does not contain a data row"
    exit 1
fi

# ============================================================
# SUMMARY
# ============================================================

log "Coverage result:"
cat "$OUT" | tee -a "$LOGFILE"

log "Finished samtools coverage for $SAMPLE"

;;

# ============================================================
# SAMTOOLS COVERAGE SUMMARY PER SAMPLE
# ============================================================
coverage_samtools_summary)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

COV_DIR="$OUTPUT_DIR/coverage_samtools"
OUT="$OUTPUT_DIR/coverage_samtools/samtools_per_sample.tsv"

log "Summarizing samtools coverage"

echo -e "sample\tmean_coverage\tbreadth" > "$OUT"

for f in "$COV_DIR"/*.samtools.tsv; do
    tail -n +2 "$f" >> "$OUT"
done

log "Created: $OUT"

;;

# ============================================================
### COVERAGE (JGI DEPTH — ARRAY OR SINGLE-SAMPLE)
# ============================================================
coverage_metabat)

set -euo pipefail
shopt -s nullglob

log "Starting JGI depth calculation (MetaBAT-compatible)"

THREADS="${THREADS:-8}"

MAP_DIR="$OUTPUT_DIR/mapping"
COV_DIR="$OUTPUT_DIR/coverage_metabat"
LOG_DIR="$OUTPUT_DIR/logs/coverage_metabat"

mkdir -p "$COV_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "MAP_DIR: $MAP_DIR"
log "COV_DIR: $COV_DIR"

# ============================================================
# TOOL CHECK
# ============================================================

command -v jgi_summarize_bam_contig_depths >/dev/null 2>&1 || {
    log "ERROR: jgi_summarize_bam_contig_depths not found"
    log "Make sure MetaBAT2 is loaded/activated in the SLURM job"
    exit 1
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    BAM="$MAP_DIR/${SAMPLE}.bam"

    [[ -s "$BAM" ]] || {
        log "ERROR: BAM not found:"
        log "$BAM"
        exit 1
    }

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    mapfile -t BAM_FILES < <(
        find "$MAP_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.bam" \
            -size +0c \
            | sort
    )

    N=${#BAM_FILES[@]}

    log "Total BAM files detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No BAM files found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    BAM="${BAM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$BAM" .bam)

fi

# ============================================================
# SAMPLE FILES
# ============================================================

DEPTH="$COV_DIR/${SAMPLE}.depth.txt"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "BAM: $BAM"
log "Depth output: $DEPTH"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# INPUT CHECK
# ============================================================

[[ -s "$BAM" ]] || {
    log "ERROR: Missing or empty BAM:"
    log "$BAM"
    exit 1
}

# ============================================================
# OPTIONAL BAM INTEGRITY CHECK
# ============================================================

if command -v samtools >/dev/null 2>&1; then

    samtools quickcheck "$BAM" || {
        log "ERROR: BAM failed samtools quickcheck:"
        log "$BAM"
        exit 1
    }

else

    log "WARNING: samtools not available -> skipping BAM integrity check"

fi

# ============================================================
# SKIP IF DONE
# ============================================================

if [[ -s "$DEPTH" && "${FORCE:-0}" != "1" ]]; then
    log "Depth already completed -> skipping"
    log "$DEPTH"
    exit 0
fi

# ============================================================
# CLEAN PARTIAL OUTPUT
# ============================================================

rm -f "$DEPTH"

# ============================================================
# RUN JGI DEPTH
# ============================================================

log "Running jgi_summarize_bam_contig_depths"

jgi_summarize_bam_contig_depths \
    --outputDepth "$DEPTH" \
    "$BAM" \
    &>> "$LOGFILE"

# ============================================================
# VALIDATE
# ============================================================

[[ -s "$DEPTH" ]] || {
    log "ERROR: Depth calculation failed"
    log "Check log:"
    log "$LOGFILE"
    exit 1
}

NLINES=$(wc -l < "$DEPTH")

if (( NLINES <= 1 )); then
    log "WARNING: Depth file contains header only"
fi

# ============================================================
# SUMMARY
# ============================================================

log "Depth preview:"
head -n 5 "$DEPTH" | tee -a "$LOGFILE"

log "Number of lines: $NLINES"
log "Finished JGI depth calculation for $SAMPLE"

;;


# ============================================================
### FILTER METABAT DEPTH
### REMOVE UNREALISTIC COVERAGE VALUES
### ARRAY OR SINGLE-SAMPLE
# ============================================================
filter_depth_metabat)

set -euo pipefail
shopt -s nullglob

log "Starting MetaBAT depth filtering"

MAX_DEPTH="${MAX_DEPTH:-10000}"

COV_DIR="$OUTPUT_DIR/coverage_metabat"
OUT_DIR="$OUTPUT_DIR/coverage_metabat_filtered"
LOG_DIR="$OUTPUT_DIR/logs/coverage_metabat_filter"

mkdir -p "$OUT_DIR" "$LOG_DIR"

log "Maximum allowed depth: $MAX_DEPTH"

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    log "Running in single-sample mode"

else

    log "Running in array mode"

    mapfile -t SAMPLE_DIRS < <(
        find "$OUTPUT_DIR/assembly" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    log "Total samples: $N"

    [[ "$N" -eq 32 ]] || {
        log "ERROR: Expected 32 sample directories but found $N"
        exit 1
    }

    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    }

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID out of range (1-$N)"
        exit 1
    fi

    SAMPLE=$(basename "${SAMPLE_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}")

fi

# ============================================================
# FILES
# ============================================================

DEPTH="$COV_DIR/${SAMPLE}.depth.txt"

OUT="$OUT_DIR/${SAMPLE}.depth.txt"

LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# ============================================================
# INPUT CHECK
# ============================================================

[[ -s "$DEPTH" ]] || {
    log "ERROR: Missing depth file"
    log "$DEPTH"
    exit 1
}

# ============================================================
# SKIP
# ============================================================

if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
    log "Filtered depth already exists -> skipping"
    exit 0
fi

# ============================================================
# FILTER
# ============================================================

TOTAL_BEFORE=$(($(wc -l < "$DEPTH") - 1))

awk -v max="$MAX_DEPTH" '
BEGIN{
    OFS="\t"
}

NR==1{
    print
    next
}

$3>=0 && $3<=max{
    print
}
' "$DEPTH" > "$OUT"

[[ -s "$OUT" ]] || {
    log "ERROR: Filtering failed"
    exit 1
}

TOTAL_AFTER=$(($(wc -l < "$OUT") - 1))

REMOVED=$((TOTAL_BEFORE-TOTAL_AFTER))

# ============================================================
# SUMMARY
# ============================================================

{
    echo "Sample: $SAMPLE"
    echo "Threshold: $MAX_DEPTH"
    echo "Contigs before: $TOTAL_BEFORE"
    echo "Contigs after:  $TOTAL_AFTER"
    echo "Removed:        $REMOVED"
} > "$LOGFILE"

log "Before filtering : $TOTAL_BEFORE"
log "After filtering  : $TOTAL_AFTER"
log "Removed          : $REMOVED"

log "Preview:"
head -5 "$OUT"

log "Finished filtering depth for $SAMPLE"

;;

# ============================================================
### CLEAN CONTIG HEADERS (MEGAHIT → BINNING-COMPATIBLE)
# ============================================================
clean_contig_headers)

set -euo pipefail
shopt -s nullglob

log "Starting contig header cleaning"

THREADS="${THREADS:-1}"

ASM_DIR="$OUTPUT_DIR/assembly"
OUT_DIR="$OUTPUT_DIR/assembly_clean"
LOG_DIR="$OUTPUT_DIR/logs/assembly_clean"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# -------------------- COLLECT SAMPLES --------------------
mapfile -t ASM_FILES < <(
    for d in "$ASM_DIR"/*; do
        [[ -d "$d" ]] || continue
        s=$(basename "$d")
        f="$d/${s}.contigs.fa"
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No assemblies found"
    exit 1
}

# MODE SELECTION

# -------------------- CASE 1: MANUAL SAMPLE --------------------
if [[ -n "${TARGET_SAMPLE:-}" ]]; then
    log "Running in MANUAL mode (TARGET_SAMPLE=$TARGET_SAMPLE)"

    ASM="$ASM_DIR/$TARGET_SAMPLE/${TARGET_SAMPLE}.contigs.fa"

    [[ -s "$ASM" ]] || {
        log "ERROR: Sample not found: $TARGET_SAMPLE"
        exit 1
    }

    SAMPLE="$TARGET_SAMPLE"

    process_sample=true

# -------------------- CASE 2: ARRAY MODE --------------------
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Running in ARRAY mode (task $SLURM_ARRAY_TASK_ID/$N)"

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID out of range (1–$N)"
        exit 1
    fi

    ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$ASM" | sed 's/.contigs.fa//')

    process_sample=true

# -------------------- CASE 3: NON-ARRAY (LOOP ALL) --------------------
else
    log "Running in NON-ARRAY mode (processing all samples)"

    process_sample=false
fi

# FUNCTION

run_clean() {
    local ASM="$1"
    local SAMPLE="$2"

    OUTDIR="$OUT_DIR/$SAMPLE"
    OUT="$OUTDIR/${SAMPLE}.contigs.clean.fa"

    log "----------------------------------------"
    log "Processing: $SAMPLE"
    log "Input:  $ASM"
    log "Output: $OUT"

    # -------------------- SKIP --------------------
    if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
        log "Clean contigs exist → skipping"
        return
    fi

    # -------------------- CLEAN --------------------
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"

    # -------------------- CLEAN HEADERS --------------------
    log "Cleaning FASTA headers"

    awk '
    /^>/ {
        split($0, a, " ")
        print a[1]
        next
    }
    {print}
    ' "$ASM" > "$OUT"

    # -------------------- VALIDATE --------------------
    [[ -s "$OUT" ]] || {
        log "ERROR: Header cleaning failed for $SAMPLE"
        exit 1
    }

    log "Finished: $SAMPLE"
}

# EXECUTION

if [[ "$process_sample" == true ]]; then
    run_clean "$ASM" "$SAMPLE"
else
    for ASM in "${ASM_FILES[@]}"; do
        SAMPLE=$(basename "$ASM" | sed 's/.contigs.fa//')
        run_clean "$ASM" "$SAMPLE"
    done
fi

log "Header cleaning step completed"
;;

# ============================================================
### BINNING (METABAT2 — CLEAN ASSEMBLY, FILTERED DEPTH)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
binning_metabat_prok)

set -euo pipefail
shopt -s nullglob

log "Starting MetaBAT2 binning (clean assembly)"

# ============================================================
# CONFIG
# ============================================================

THREADS="${THREADS:-8}"
MIN_LEN="${MIN_LEN:-1500}"

ASM_DIR="$OUTPUT_DIR/assembly_clean"

COV_DIR="$OUTPUT_DIR/coverage_metabat"
COV_FILT_DIR="$OUTPUT_DIR/coverage_metabat_filtered"

BIN_DIR="$OUTPUT_DIR/binning_metabat_prok"
LOG_DIR="$OUTPUT_DIR/logs/binning_metabat_prok"

mkdir -p "$BIN_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "MIN_LEN: $MIN_LEN"
log "ASM_DIR: $ASM_DIR"
log "COV_DIR: $COV_DIR"
log "COV_FILT_DIR: $COV_FILT_DIR"
log "BIN_DIR: $BIN_DIR"

# ============================================================
# INPUT DIRECTORY CHECK
# ============================================================

[[ -d "$ASM_DIR" ]] || {
    log "ERROR: Clean assembly directory not found:"
    log "$ASM_DIR"
    exit 1
}

# ============================================================
# TOOL CHECKS
# ============================================================

command -v metabat2 >/dev/null 2>&1 || {
    log "ERROR: metabat2 not found"
    exit 1
}

command -v seqkit >/dev/null 2>&1 || {
    log "ERROR: seqkit not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_metabat() {

    local SAMPLE="$1"

    local CONTIGS="$ASM_DIR/$SAMPLE/${SAMPLE}.contigs.clean.fa"

    local DEPTH_FILT="$COV_FILT_DIR/${SAMPLE}.depth.txt"
    local DEPTH_RAW="$COV_DIR/${SAMPLE}.depth.txt"

    local DEPTH=""

    local OUTDIR="$BIN_DIR/$SAMPLE"

    local FILTERED_CONTIGS="$OUTDIR/${SAMPLE}.min${MIN_LEN}.fa"

    local DEPTH_LIST="$OUTDIR/depth.list"
    local CONTIG_LIST="$OUTDIR/contigs.list"
    local SHARED_LIST="$OUTDIR/shared.list"

    local MATCHED_CONTIGS="$OUTDIR/${SAMPLE}.metabat.input.fa"
    local MATCHED_DEPTH="$OUTDIR/${SAMPLE}.metabat.depth.txt"

    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local N_SHARED=0
    local N_BINS=0
    local EXIT_CODE=0

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "Clean assembly: $CONTIGS"
    log "Output directory: $OUTDIR"

    # ========================================================
    # INPUT CHECK: CONTIGS
    # ========================================================

    if [[ ! -e "$CONTIGS" ]]; then
        log "ERROR: Clean assembly does not exist:"
        log "$CONTIGS"
        return 1
    fi

    if [[ ! -s "$CONTIGS" ]]; then
        log "No clean contigs for $SAMPLE -> MetaBAT2 not required"
        return 0
    fi

    # ========================================================
    # SELECT DEPTH FILE
    # ========================================================

    if [[ -s "$DEPTH_FILT" ]]; then

        DEPTH="$DEPTH_FILT"

        log "Using FILTERED depth:"
        log "$DEPTH"

    elif [[ -s "$DEPTH_RAW" ]]; then

        DEPTH="$DEPTH_RAW"

        log "WARNING: Filtered depth not found"
        log "Using UNFILTERED depth:"
        log "$DEPTH"

    else

        log "ERROR: No depth file found for $SAMPLE"
        log "Checked:"
        log "$DEPTH_FILT"
        log "$DEPTH_RAW"

        return 1
    fi

    # ========================================================
    # SKIP IF COMPLETED
    # ========================================================

    if [[ "${FORCE:-0}" != "1" ]]; then

        if [[ -f "$OUTDIR/METABAT_COMPLETE" ]]; then
            log "MetaBAT2 already completed -> skipping"
            return 0
        fi

        if [[ -f "$OUTDIR/NO_BINS" ]]; then
            log "MetaBAT2 previously completed with no bins -> skipping"
            return 0
        fi

        if [[ -d "$OUTDIR" ]] &&
           compgen -G "$OUTDIR/bin*.fa" >/dev/null; then

            log "MetaBAT2 bins already exist -> skipping"
            return 0
        fi

    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then

        log "Removing previous or incomplete MetaBAT2 directory:"
        log "$OUTDIR"

        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    rm -f "$LOGFILE"

    # ========================================================
    # FILTER CONTIGS BY LENGTH
    # ========================================================

    log "Filtering contigs"
    log "Minimum contig length: ${MIN_LEN} bp"

    seqkit seq \
        -m "$MIN_LEN" \
        "$CONTIGS" \
        > "$FILTERED_CONTIGS"

    if [[ ! -s "$FILTERED_CONTIGS" ]]; then

        log "No contigs >= ${MIN_LEN} bp"
        log "Nothing to bin for $SAMPLE"

        touch "$OUTDIR/NO_BINS"

        return 0
    fi

    # ========================================================
    # CONTIG / DEPTH CONSISTENCY CHECK
    # ========================================================

    log "Checking contig-depth consistency"

    # Depth contig IDs
    awk -F'\t' '
        NR > 1 && NF > 0 {
            print $1
        }
    ' "$DEPTH" \
        | sort -u \
        > "$DEPTH_LIST"

    # FASTA contig IDs
    awk '
        /^>/ {
            sub(/^>/,"")
            split($0,a,/[\t ]/)
            print a[1]
        }
    ' "$FILTERED_CONTIGS" \
        | sort -u \
        > "$CONTIG_LIST"

    # Shared contigs
    comm -12 \
        "$DEPTH_LIST" \
        "$CONTIG_LIST" \
        > "$SHARED_LIST"

    N_SHARED=$(wc -l < "$SHARED_LIST")

    log "Matching contigs between FASTA and depth: $N_SHARED"

    if (( N_SHARED == 0 )); then

        log "ERROR: No overlap between assembly contigs and depth file"
        log "MetaBAT2 cannot proceed"

        return 1

    elif (( N_SHARED < 100 )); then

        log "WARNING: Low contig-depth overlap ($N_SHARED)"

    fi

    # ========================================================
    # CREATE MATCHED FASTA
    # ========================================================

    log "Creating MetaBAT2 FASTA containing only shared contigs"

    seqkit grep \
        -f "$SHARED_LIST" \
        "$FILTERED_CONTIGS" \
        > "$MATCHED_CONTIGS"

    if [[ ! -s "$MATCHED_CONTIGS" ]]; then
        log "ERROR: Matched FASTA is empty"
        return 1
    fi

    # ========================================================
    # CREATE MATCHED DEPTH FILE
    # ========================================================

    log "Creating matching MetaBAT2 depth file"

    awk -F'\t' '
        BEGIN {
            OFS="\t"
        }

        NR==FNR {
            keep[$1]=1
            next
        }

        FNR==1 {
            print
            next
        }

        ($1 in keep) {
            print
        }
    ' "$SHARED_LIST" "$DEPTH" \
        > "$MATCHED_DEPTH"

    if [[ ! -s "$MATCHED_DEPTH" ]]; then
        log "ERROR: Matched depth file is empty"
        return 1
    fi

    # ========================================================
    # RUN METABAT2
    # ========================================================

    log "Running MetaBAT2"
    log "Input FASTA: $MATCHED_CONTIGS"
    log "Input depth: $MATCHED_DEPTH"

    set +e

    metabat2 \
        -i "$MATCHED_CONTIGS" \
        -a "$MATCHED_DEPTH" \
        -o "$OUTDIR/bin" \
        -t "$THREADS" \
        -m "$MIN_LEN" \
        &>> "$LOGFILE"

    EXIT_CODE=$?

    set -e

    log "MetaBAT2 exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if (( EXIT_CODE != 0 )); then

        log "ERROR: MetaBAT2 failed for $SAMPLE"
        log "Last 50 lines of log:"

        tail -n 50 "$LOGFILE" || true

        return "$EXIT_CODE"
    fi

    # ========================================================
    # VALIDATE / COUNT BINS
    # ========================================================

    N_BINS=$(find "$OUTDIR" \
        -maxdepth 1 \
        -type f \
        -name "bin*.fa" \
        | wc -l)

    if (( N_BINS == 0 )); then

        log "MetaBAT2 completed successfully but produced no bins"

        touch "$OUTDIR/NO_BINS"

    else

        log "Bins produced: $N_BINS"

        touch "$OUTDIR/METABAT_COMPLETE"

    fi

    log "Finished MetaBAT2 for $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"
    log "Sample: $SAMPLE"

    run_metabat "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$ASM_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.clean.fa" \
            -size +0c \
            | sort
    )

    N=${#ASM_FILES[@]}

    log "Total cleaned assemblies detected: $N"

    if (( N == 0 )); then
        log "ERROR: No cleaned assemblies found"
        exit 1
    fi

    # -------------------- ARRAY CHECK --------------------

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------

    CONTIGS="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$CONTIGS" .contigs.clean.fa)

    log "Array task ${SLURM_ARRAY_TASK_ID} / $N -> $SAMPLE"

    run_metabat "$SAMPLE"

fi

log "MetaBAT2 prokaryotic binning step completed"

;;

# ============================================================
### CLEAN COVERM (MAKE MAXBIN-COMPATIBLE)
# Uses output from coverage_coverm
# ============================================================
clean_coverm)

set -euo pipefail
shopt -s nullglob

log "Starting CoverM cleaning for MaxBin"

# ============================================================
# CONFIG
# ============================================================

COV_DIR="$OUTPUT_DIR/coverage_coverm"
ASM_CLEAN_DIR="$OUTPUT_DIR/assembly_clean"

OUT_DIR="$OUTPUT_DIR/coverage_coverm_clean"
LOG_DIR="$OUTPUT_DIR/logs/coverage_coverm_clean"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Backward compatibility
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

# ============================================================
# FUNCTION
# ============================================================

run_clean_coverm() {

    local SAMPLE="$1"

    local COV="$COV_DIR/${SAMPLE}.coverage.tsv"

    local CLEAN_CONTIGS="$ASM_CLEAN_DIR/$SAMPLE/${SAMPLE}.contigs.clean.fa"

    local OUT="$OUT_DIR/${SAMPLE}.coverm.clean.tsv"

    local CONTIG_LIST="$OUT_DIR/${SAMPLE}.contigs.list"

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "Coverage: $COV"
    log "Clean contigs: $CLEAN_CONTIGS"
    log "Output: $OUT"

    # ========================================================
    # INPUT CHECKS
    # ========================================================

    [[ -s "$COV" ]] || {
        log "ERROR: Coverage file missing:"
        log "$COV"
        return 1
    }

    [[ -s "$CLEAN_CONTIGS" ]] || {
        log "ERROR: Cleaned assembly missing:"
        log "$CLEAN_CONTIGS"
        return 1
    }

    # ========================================================
    # SKIP IF DONE
    # ========================================================

    if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then

        log "Clean CoverM already exists -> skipping"

        return 0
    fi

    # ========================================================
    # EXTRACT CLEAN CONTIG IDs
    # ========================================================

    awk '
    /^>/ {
        sub(/^>/,"")
        split($0,a,/[\t ]/)
        print a[1]
    }
    ' "$CLEAN_CONTIGS" \
        | sort -u \
        > "$CONTIG_LIST"

    N_CONTIGS=$(wc -l < "$CONTIG_LIST")

    log "Clean contigs detected: $N_CONTIGS"

    if (( N_CONTIGS == 0 )); then
        log "ERROR: No cleaned contigs detected"
        return 1
    fi

    # ========================================================
    # FILTER COVERAGE
    #
    # coverage_coverm output:
    #
    # $1 = contig
    # $2 = mean_coverage
    # $3 = covered_fraction
    # $4 = length
    # $5 = coverage_weighted
    #
    # MaxBin requires:
    #
    # contig_ID    abundance
    #
    # Therefore use mean_coverage ($2)
    # ========================================================

    awk '
    BEGIN {
        FS=OFS="\t"
    }

    NR==FNR {
        keep[$1]=1
        next
    }

    FNR==1 {
        next
    }

    ($1 in keep) && $2 != "NA" && ($2+0) > 0 {
        print $1, $2
    }
    ' "$CONTIG_LIST" "$COV" > "$OUT"

    # ========================================================
    # VALIDATE
    # ========================================================

    if [[ ! -s "$OUT" ]]; then

        log "ERROR: No matching coverage entries produced"

        return 1
    fi

    N_SHARED=$(wc -l < "$OUT")

    log "Contigs with coverage: $N_SHARED"

    if (( N_SHARED < 10 )); then
        log "WARNING: Only $N_SHARED contigs have usable coverage"
    fi

    # ========================================================
    # PREVIEW
    # ========================================================

    log "Preview:"
    head -n 5 "$OUT"

    log "Finished clean CoverM for $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"

    run_clean_coverm "$SAMPLE"

elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t COV_FILES < <(
        find "$COV_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.coverage.tsv" \
            | sort
    )

    N=${#COV_FILES[@]}

    log "Coverage files detected: $N"

    if (( N == 0 )); then
        log "ERROR: No coverage files found"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then

        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"

        exit 1
    fi

    COV="${COV_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$COV" .coverage.tsv)

    log "Array task ${SLURM_ARRAY_TASK_ID} -> $SAMPLE"

    run_clean_coverm "$SAMPLE"

else

    # --------------------------------------------------------
    # LOCAL SEQUENTIAL MODE
    # --------------------------------------------------------

    log "Running in local sequential mode"

    mapfile -t COV_FILES < <(
        find "$COV_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.coverage.tsv" \
            | sort
    )

    N=${#COV_FILES[@]}

    if (( N == 0 )); then
        log "ERROR: No coverage files found"
        exit 1
    fi

    for COV in "${COV_FILES[@]}"; do

        SAMPLE=$(basename "$COV" .coverage.tsv)

        run_clean_coverm "$SAMPLE"

    done

fi

log "CoverM cleaning completed"

;;

# ============================================================
### BINNING (MAXBIN2 — CLEAN ASSEMBLY + CLEAN COVERM)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
binning_maxbin_prok)

set -euo pipefail
shopt -s nullglob

log "Starting MaxBin2 (clean assembly + clean CoverM)"

# ============================================================
# CONFIG
# ============================================================

THREADS="${THREADS:-8}"
MIN_LEN="${MIN_LEN:-1500}"

ASM_DIR="$OUTPUT_DIR/assembly_clean"
COV_DIR="$OUTPUT_DIR/coverage_coverm_clean"

BIN_DIR="$OUTPUT_DIR/binning_maxbin_prok"
LOG_DIR="$OUTPUT_DIR/logs/binning_maxbin_prok"

mkdir -p "$BIN_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "MIN_LEN: $MIN_LEN"
log "ASM_DIR: $ASM_DIR"
log "COV_DIR: $COV_DIR"
log "BIN_DIR: $BIN_DIR"

# ============================================================
# INPUT DIRECTORY CHECKS
# ============================================================

[[ -d "$ASM_DIR" ]] || {
    log "ERROR: Clean assembly directory not found:"
    log "$ASM_DIR"
    exit 1
}

[[ -d "$COV_DIR" ]] || {
    log "ERROR: Clean CoverM directory not found:"
    log "$COV_DIR"
    exit 1
}

# ============================================================
# TOOL CHECKS
# ============================================================

command -v run_MaxBin.pl >/dev/null 2>&1 || {
    log "ERROR: run_MaxBin.pl not found"
    exit 1
}

command -v seqkit >/dev/null 2>&1 || {
    log "ERROR: seqkit not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_maxbin() {

    local SAMPLE="$1"

    local CONTIGS="$ASM_DIR/$SAMPLE/${SAMPLE}.contigs.clean.fa"

    local COVERM="$COV_DIR/${SAMPLE}.coverm.clean.tsv"

    local OUTDIR="$BIN_DIR/$SAMPLE"

    local FILTERED_CONTIGS="$OUTDIR/${SAMPLE}.min${MIN_LEN}.fa"
    local CONTIG_LIST="$OUTDIR/contigs.list"
    local ABUND_LIST="$OUTDIR/abund.list"
    local SHARED_LIST="$OUTDIR/shared.list"

    local MATCHED_CONTIGS="$OUTDIR/${SAMPLE}.maxbin.input.fa"
    local ABUND="$OUTDIR/${SAMPLE}.maxbin.tsv"

    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local N_SHARED=0
    local N_BINS=0
    local EXIT_CODE=0

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "Clean assembly: $CONTIGS"
    log "Clean CoverM: $COVERM"
    log "Output directory: $OUTDIR"

    # ========================================================
    # INPUT CHECKS
    # ========================================================

    if [[ ! -e "$CONTIGS" ]]; then
        log "ERROR: Clean assembly does not exist:"
        log "$CONTIGS"
        return 1
    fi

    if [[ ! -s "$CONTIGS" ]]; then
        log "No clean contigs for $SAMPLE -> MaxBin2 not required"
        return 0
    fi

    if [[ ! -e "$COVERM" ]]; then
        log "ERROR: Clean CoverM file does not exist:"
        log "$COVERM"
        return 1
    fi

    if [[ ! -s "$COVERM" ]]; then
        log "No usable coverage values for $SAMPLE -> MaxBin2 not required"
        return 0
    fi

    # ========================================================
    # SKIP IF COMPLETED
    # ========================================================

    if [[ "${FORCE:-0}" != "1" ]]; then

        if [[ -f "$OUTDIR/MAXBIN_COMPLETE" ]]; then
            log "MaxBin2 already completed -> skipping"
            return 0
        fi

        if [[ -f "$OUTDIR/NO_BINS" ]]; then
            log "MaxBin2 previously completed with no bins -> skipping"
            return 0
        fi

        if [[ -d "$OUTDIR" ]] &&
           compgen -G "$OUTDIR/bin.*.fasta" >/dev/null; then

            log "MaxBin2 bins already exist -> skipping"
            return 0
        fi

    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then

        log "Removing previous or incomplete MaxBin2 output:"
        log "$OUTDIR"

        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    rm -f "$LOGFILE"

    # ========================================================
    # FILTER CONTIGS BY LENGTH
    # ========================================================

    log "Filtering contigs"
    log "Minimum contig length: ${MIN_LEN} bp"

    seqkit seq \
        -m "$MIN_LEN" \
        "$CONTIGS" \
        > "$FILTERED_CONTIGS"

    if [[ ! -s "$FILTERED_CONTIGS" ]]; then

        log "No contigs >= ${MIN_LEN} bp"
        log "Nothing to bin for $SAMPLE"

        touch "$OUTDIR/NO_BINS"

        return 0
    fi

    # ========================================================
    # CHECK CONTIG / ABUNDANCE OVERLAP
    # ========================================================

    log "Checking contig-abundance consistency"

    # FASTA IDs: first token only
    awk '
        /^>/ {
            sub(/^>/,"")
            split($0,a,/[\t ]/)
            print a[1]
        }
    ' "$FILTERED_CONTIGS" \
        | sort -u \
        > "$CONTIG_LIST"

    # CoverM abundance IDs
    cut -f1 "$COVERM" \
        | sort -u \
        > "$ABUND_LIST"

    comm -12 \
        "$CONTIG_LIST" \
        "$ABUND_LIST" \
        > "$SHARED_LIST"

    N_SHARED=$(wc -l < "$SHARED_LIST")

    log "Matching contigs: $N_SHARED"

    if (( N_SHARED == 0 )); then

        log "ERROR: No overlap between cleaned assembly and abundance file"

        return 1

    elif (( N_SHARED < 100 )); then

        log "WARNING: Low contig-abundance overlap ($N_SHARED)"

    fi

    # ========================================================
    # CREATE MATCHED FASTA
    # ========================================================

    log "Creating MaxBin FASTA containing shared contigs"

    seqkit grep \
        -f "$SHARED_LIST" \
        "$FILTERED_CONTIGS" \
        > "$MATCHED_CONTIGS"

    if [[ ! -s "$MATCHED_CONTIGS" ]]; then
        log "ERROR: Matched MaxBin FASTA is empty"
        return 1
    fi

    # ========================================================
    # CREATE MATCHED ABUNDANCE FILE
    # ========================================================

    log "Preparing abundance file for MaxBin"

    awk '
        BEGIN {
            FS=OFS="\t"
        }

        NR==FNR {
            keep[$1]=1
            next
        }

        ($1 in keep) {
            print $1, $2
        }
    ' "$SHARED_LIST" "$COVERM" \
        > "$ABUND"

    if [[ ! -s "$ABUND" ]]; then
        log "ERROR: MaxBin abundance file is empty"
        return 1
    fi

    log "Abundance entries: $(wc -l < "$ABUND")"

    # ========================================================
    # RUN MAXBIN2
    # ========================================================

    log "Running MaxBin2"
    log "Input FASTA: $MATCHED_CONTIGS"
    log "Abundance: $ABUND"

    set +e

    run_MaxBin.pl \
        -contig "$MATCHED_CONTIGS" \
        -abund "$ABUND" \
        -out "$OUTDIR/bin" \
        -thread "$THREADS" \
        &>> "$LOGFILE"

    EXIT_CODE=$?

    set -e

    log "MaxBin2 exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if (( EXIT_CODE != 0 )); then

        log "ERROR: MaxBin2 failed for $SAMPLE"
        log "Last 50 lines of MaxBin log:"

        tail -n 50 "$LOGFILE" || true

        return "$EXIT_CODE"
    fi

    # ========================================================
    # VALIDATE / COUNT BINS
    # ========================================================

    N_BINS=$(find "$OUTDIR" \
        -maxdepth 1 \
        -type f \
        -name "bin.*.fasta" \
        | wc -l)

    if (( N_BINS == 0 )); then

        log "MaxBin2 completed successfully but produced no bins"

        touch "$OUTDIR/NO_BINS"

    else

        log "Bins produced: $N_BINS"

        touch "$OUTDIR/MAXBIN_COMPLETE"

    fi

    log "Finished MaxBin2 for $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility with older SLURM scripts
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"
    log "Sample: $SAMPLE"

    run_maxbin "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$ASM_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.clean.fa" \
            -size +0c \
            | sort
    )

    N=${#ASM_FILES[@]}

    log "Total cleaned assemblies detected: $N"

    if (( N == 0 )); then
        log "ERROR: No cleaned assemblies found"
        exit 1
    fi

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    CONTIGS="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$CONTIGS" .contigs.clean.fa)

    log "Array task ${SLURM_ARRAY_TASK_ID} / $N -> $SAMPLE"

    run_maxbin "$SAMPLE"

fi

log "MaxBin2 prokaryotic binning step completed"

;;

# ============================================================
### BIN REFINEMENT (DAS TOOL — PROK)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
das_tool_prok)

set -euo pipefail
shopt -s nullglob

log "Starting DAS Tool (prok bin refinement)"

# ============================================================
# CONFIG
# ============================================================

THREADS="${THREADS:-16}"
MIN_LEN="${MIN_LEN:-1500}"

MB_DIR="$OUTPUT_DIR/binning_metabat_prok"
MX_DIR="$OUTPUT_DIR/binning_maxbin_prok"
ASM_CLEAN_DIR="$OUTPUT_DIR/assembly_clean"

DAS_DIR="$OUTPUT_DIR/das_tool_prok"
LOG_DIR="$OUTPUT_DIR/logs/das_tool_prok"

mkdir -p "$DAS_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "MIN_LEN: $MIN_LEN"
log "MetaBAT directory: $MB_DIR"
log "MaxBin directory: $MX_DIR"
log "Clean assembly directory: $ASM_CLEAN_DIR"
log "DAS Tool directory: $DAS_DIR"

# ============================================================
# INPUT DIRECTORY CHECKS
# ============================================================

[[ -d "$ASM_CLEAN_DIR" ]] || {
    log "ERROR: Clean assembly directory not found:"
    log "$ASM_CLEAN_DIR"
    exit 1
}

[[ -d "$MB_DIR" ]] || {
    log "ERROR: MetaBAT directory not found:"
    log "$MB_DIR"
    exit 1
}

[[ -d "$MX_DIR" ]] || {
    log "ERROR: MaxBin directory not found:"
    log "$MX_DIR"
    exit 1
}

# ============================================================
# TOOL CHECKS
# ============================================================

command -v DAS_Tool >/dev/null 2>&1 || {
    log "ERROR: DAS_Tool not found"
    exit 1
}

command -v seqkit >/dev/null 2>&1 || {
    log "ERROR: seqkit not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_dastool() {

    local SAMPLE="$1"

    local MB="$MB_DIR/$SAMPLE"
    local MX="$MX_DIR/$SAMPLE"

    local ASM_CLEAN="$ASM_CLEAN_DIR/$SAMPLE/${SAMPLE}.contigs.clean.fa"

    local OUTDIR="$DAS_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}"

    local FILTERED_CONTIGS="$OUTDIR/${SAMPLE}.min${MIN_LEN}.fa"
    local CONTIG_LIST="$OUTDIR/contigs.list"

    local MB_MAP_RAW="$OUTDIR/metabat_scaffolds2bin.raw.tsv"
    local MX_MAP_RAW="$OUTDIR/maxbin_scaffolds2bin.raw.tsv"

    local MB_MAP="$OUTDIR/metabat_scaffolds2bin.tsv"
    local MX_MAP="$OUTDIR/maxbin_scaffolds2bin.tsv"

    local SUMMARY_FILE="${PREFIX}_DASTool_summary.tsv"
    local CONTIG2BIN_FILE="${PREFIX}_DASTool_contig2bin.tsv"
    local BINS_DIR="${PREFIX}_DASTool_bins"

    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local INPUT_FILES=""
    local LABELS=""
    local N_INPUTS=0
    local N_BINS=0
    local N_SHARED=0
    local EXIT_CODE=0

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "Assembly: $ASM_CLEAN"
    log "MetaBAT: $MB"
    log "MaxBin: $MX"
    log "Output: $OUTDIR"

    # ========================================================
    # INPUT CHECK: ASSEMBLY
    # ========================================================

    if [[ ! -e "$ASM_CLEAN" ]]; then
        log "ERROR: Clean assembly does not exist:"
        log "$ASM_CLEAN"
        return 1
    fi

    if [[ ! -s "$ASM_CLEAN" ]]; then
        log "No clean contigs for $SAMPLE -> DAS Tool not required"
        return 0
    fi

    # ========================================================
    # SKIP IF ALREADY COMPLETE
    # ========================================================

    if [[ "${FORCE:-0}" != "1" ]]; then

        if [[ -f "$OUTDIR/DASTOOL_COMPLETE" ]]; then
            log "DAS Tool already completed -> skipping"
            return 0
        fi

        if [[ -f "$OUTDIR/NO_BINS" ]]; then
            log "DAS Tool previously completed with no bins -> skipping"
            return 0
        fi

        # Backward compatibility with older completed runs
        if [[ -s "$SUMMARY_FILE" &&
              -s "$CONTIG2BIN_FILE" &&
              -d "$BINS_DIR" ]]; then

            N_BINS=$(find "$BINS_DIR" \
                -maxdepth 1 \
                -type f \
                -name "*.fa" \
                | wc -l)

            if (( N_BINS > 0 )); then
                log "Existing DAS Tool result detected ($N_BINS bins) -> skipping"
                return 0
            fi
        fi

    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then

        log "Removing previous or incomplete DAS Tool output:"
        log "$OUTDIR"

        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    rm -f "$LOGFILE"

    # ========================================================
    # FILTER ASSEMBLY TO SAME MINIMUM LENGTH AS BINNING
    # ========================================================

    log "Filtering contigs"
    log "Minimum length: ${MIN_LEN} bp"

    seqkit seq \
        -m "$MIN_LEN" \
        "$ASM_CLEAN" \
        > "$FILTERED_CONTIGS"

    if [[ ! -s "$FILTERED_CONTIGS" ]]; then

        log "No contigs >= ${MIN_LEN} bp"
        log "Nothing available for DAS Tool"

        touch "$OUTDIR/NO_BINS"

        return 0
    fi

    # ========================================================
    # ASSEMBLY CONTIG IDS
    # ========================================================

    awk '
        /^>/ {
            sub(/^>/,"")
            split($0,a,/[\t ]/)
            print a[1]
        }
    ' "$FILTERED_CONTIGS" \
        | sort -u \
        > "$CONTIG_LIST"

    # ========================================================
    # CREATE METABAT SCAFFOLD2BIN
    # ========================================================

    : > "$MB_MAP_RAW"

    if [[ -d "$MB" ]]; then

        for f in "$MB"/bin*.fa; do

            [[ -s "$f" ]] || continue

            BIN=$(basename "$f" .fa)

            awk -v bin="$BIN" '
                /^>/ {
                    sub(/^>/,"")
                    split($0,a,/[\t ]/)
                    print a[1] "\t" bin
                }
            ' "$f"

        done \
            | sort -u \
            > "$MB_MAP_RAW"

    fi

    # ========================================================
    # CREATE MAXBIN SCAFFOLD2BIN
    # ========================================================

    : > "$MX_MAP_RAW"

    if [[ -d "$MX" ]]; then

        for f in "$MX"/bin.*.fasta; do

            [[ -s "$f" ]] || continue

            BIN=$(basename "$f" .fasta)

            awk -v bin="$BIN" '
                /^>/ {
                    sub(/^>/,"")
                    split($0,a,/[\t ]/)
                    print a[1] "\t" bin
                }
            ' "$f"

        done \
            | sort -u \
            > "$MX_MAP_RAW"

    fi

    # ========================================================
    # FILTER MAPS TO CONTIGS PRESENT IN DAS TOOL ASSEMBLY
    # ========================================================

    awk '
        BEGIN {FS=OFS="\t"}

        NR==FNR {
            keep[$1]=1
            next
        }

        ($1 in keep) {
            print
        }
    ' "$CONTIG_LIST" "$MB_MAP_RAW" \
        > "$MB_MAP"

    awk '
        BEGIN {FS=OFS="\t"}

        NR==FNR {
            keep[$1]=1
            next
        }

        ($1 in keep) {
            print
        }
    ' "$CONTIG_LIST" "$MX_MAP_RAW" \
        > "$MX_MAP"

    # ========================================================
    # BUILD DAS TOOL INPUT LIST
    # ========================================================

    if [[ -s "$MB_MAP" ]]; then

        INPUT_FILES="$MB_MAP"
        LABELS="metabat"
        N_INPUTS=$((N_INPUTS + 1))

        log "MetaBAT scaffold2bin entries: $(wc -l < "$MB_MAP")"
    else
        log "WARNING: No usable MetaBAT bins for $SAMPLE"
    fi

    if [[ -s "$MX_MAP" ]]; then

        if [[ -n "$INPUT_FILES" ]]; then
            INPUT_FILES="${INPUT_FILES},${MX_MAP}"
            LABELS="${LABELS},maxbin"
        else
            INPUT_FILES="$MX_MAP"
            LABELS="maxbin"
        fi

        N_INPUTS=$((N_INPUTS + 1))

        log "MaxBin scaffold2bin entries: $(wc -l < "$MX_MAP")"
    else
        log "WARNING: No usable MaxBin bins for $SAMPLE"
    fi

    # ========================================================
    # NO BINS FROM EITHER METHOD
    # ========================================================

    if (( N_INPUTS == 0 )); then

        log "Neither MetaBAT nor MaxBin produced usable bins"
        log "DAS Tool not required"

        touch "$OUTDIR/NO_BINS"

        return 0
    fi

    if (( N_INPUTS == 1 )); then
        log "WARNING: Only one binning method produced usable bins"
        log "DAS Tool will run with the available bin set only"
    fi

    # ========================================================
    # CONSISTENCY CHECK
    # ========================================================

    log "Checking contig consistency"

    cat "$MB_MAP" "$MX_MAP" \
        | awk -F'\t' 'NF > 0 {print $1}' \
        | sort -u \
        > "$OUTDIR/bin_contigs.list"

    comm -12 \
        "$CONTIG_LIST" \
        "$OUTDIR/bin_contigs.list" \
        > "$OUTDIR/shared.list"

    N_SHARED=$(wc -l < "$OUTDIR/shared.list")

    log "Shared contigs between bins and assembly: $N_SHARED"

    if (( N_SHARED == 0 )); then

        log "ERROR: No overlap between binned contigs and assembly"

        return 1

    elif (( N_SHARED < 100 )); then

        log "WARNING: Low overlap between bins and assembly ($N_SHARED)"

    fi

    # ========================================================
    # RUN DAS TOOL
    # ========================================================

    log "Running DAS Tool"
    log "Input bin sets: $INPUT_FILES"
    log "Labels: $LABELS"

    set +e

    DAS_Tool \
        -i "$INPUT_FILES" \
        -c "$FILTERED_CONTIGS" \
        -o "$PREFIX" \
        --labels "$LABELS" \
        --threads "$THREADS" \
        --search_engine diamond \
        --write_bins \
        --score_threshold 0.5 \
        &>> "$LOGFILE"

    EXIT_CODE=$?

    set -e

    log "DAS Tool exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if (( EXIT_CODE != 0 )); then

        log "ERROR: DAS Tool failed for $SAMPLE"
        log "Last 50 lines of DAS Tool log:"

        tail -n 50 "$LOGFILE" || true

        return "$EXIT_CODE"
    fi

    # ========================================================
    # VALIDATION
    # ========================================================

    if [[ -d "$BINS_DIR" ]]; then

        N_BINS=$(find "$BINS_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.fa" \
            | wc -l)

    else

        N_BINS=0
    fi

    if (( N_BINS > 0 )); then

        if [[ ! -s "$SUMMARY_FILE" ]]; then
            log "ERROR: DAS Tool bins exist but summary file is missing"
            return 1
        fi

        if [[ ! -s "$CONTIG2BIN_FILE" ]]; then
            log "ERROR: DAS Tool bins exist but contig2bin file is missing"
            return 1
        fi

        log "SUCCESS: DAS Tool produced $N_BINS refined bins"

        touch "$OUTDIR/DASTOOL_COMPLETE"

    else

        log "DAS Tool completed successfully but selected no refined bins"
        log "This is treated as a valid zero-bin result"

        touch "$OUTDIR/NO_BINS"

    fi

    log "Finished DAS Tool for $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility with older scripts
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"
    log "Sample: $SAMPLE"

    run_dastool "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # Stable list based on cleaned assemblies.
    # Do NOT build an array from only incomplete samples.
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$ASM_CLEAN_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.clean.fa" \
            -size +0c \
            | sort
    )

    N=${#ASM_FILES[@]}

    log "Total cleaned assemblies detected: $N"

    if (( N == 0 )); then
        log "ERROR: No cleaned assemblies found"
        exit 1
    fi

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$ASM" .contigs.clean.fa)

    log "Array task ${SLURM_ARRAY_TASK_ID} / $N -> $SAMPLE"

    run_dastool "$SAMPLE"

fi

log "DAS Tool prokaryotic refinement step completed"

;;

# ============================================================
### CHECKM (PROK QUALITY — FROM DAS TOOL)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
checkm_prok)

set -euo pipefail
shopt -s nullglob

log "Starting CheckM (from DAS Tool refined bins)"

# ============================================================
# CONFIG
# ============================================================

THREADS="${THREADS:-16}"
MIN_SIZE="${MIN_SIZE:-50000}"

DAS_DIR="$OUTPUT_DIR/das_tool_prok"
ASM_CLEAN_DIR="$OUTPUT_DIR/assembly_clean"

CHECKM_DIR="$OUTPUT_DIR/checkm_prok"
LOG_DIR="$OUTPUT_DIR/logs/checkm_prok"

mkdir -p "$CHECKM_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "MIN_SIZE: $MIN_SIZE"
log "DAS_DIR: $DAS_DIR"
log "CHECKM_DIR: $CHECKM_DIR"

# ============================================================
# INPUT DIRECTORY CHECKS
# ============================================================

[[ -d "$DAS_DIR" ]] || {
    log "ERROR: DAS Tool directory not found:"
    log "$DAS_DIR"
    exit 1
}

[[ -d "$ASM_CLEAN_DIR" ]] || {
    log "ERROR: Clean assembly directory not found:"
    log "$ASM_CLEAN_DIR"
    exit 1
}

# ============================================================
# TOOL CHECK
# ============================================================

command -v checkm >/dev/null 2>&1 || {
    log "ERROR: CheckM not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_checkm() {

    local SAMPLE="$1"

    local DAS_SAMPLE_DIR="$DAS_DIR/$SAMPLE"
    local BIN_DIR="$DAS_SAMPLE_DIR/${SAMPLE}_DASTool_bins"

    local OUTDIR="$CHECKM_DIR/$SAMPLE"
    local INPUT_DIR="$OUTDIR/input"

    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local QA_FILE="$OUTDIR/${SAMPLE}_qa.tsv"
    local LINEAGE_FILE="$OUTDIR/lineage.ms"
    local BIN_STATS="$OUTDIR/storage/bin_stats_ext.tsv"

    local N_KEPT=0
    local EXIT_CODE=0
    local QA_EXIT=0

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "DAS Tool sample directory: $DAS_SAMPLE_DIR"
    log "DAS Tool bins: $BIN_DIR"
    log "CheckM output: $OUTDIR"

    # ========================================================
    # CHECK UPSTREAM DAS TOOL STATUS
    # ========================================================

    if [[ -f "$DAS_SAMPLE_DIR/NO_BINS" ]]; then

        log "DAS Tool produced no refined bins for $SAMPLE"
        log "CheckM not required"

        mkdir -p "$OUTDIR"
        touch "$OUTDIR/NO_BINS"

        return 0
    fi

    # ========================================================
    # CHECK FOR COMPLETED CHECKM OUTPUT
    # ========================================================

    if [[ "${FORCE:-0}" != "1" ]]; then

        if [[ -f "$OUTDIR/CHECKM_COMPLETE" ]]; then
            log "CheckM already completed -> skipping"
            return 0
        fi

        if [[ -f "$OUTDIR/NO_BINS" ]]; then
            log "No bins previously available for CheckM -> skipping"
            return 0
        fi

        # Backward compatibility with completed runs
        if [[ -s "$BIN_STATS" &&
              -s "$LINEAGE_FILE" &&
              -s "$QA_FILE" ]]; then

            log "Existing complete CheckM outputs detected -> skipping"
            return 0
        fi

    fi

    # ========================================================
    # CHECK DAS TOOL BIN DIRECTORY
    # ========================================================

    if [[ ! -d "$BIN_DIR" ]]; then

        log "ERROR: DAS Tool bins directory missing:"
        log "$BIN_DIR"

        log "If DAS Tool legitimately produced no bins, a NO_BINS marker"
        log "should exist in:"
        log "$DAS_SAMPLE_DIR"

        return 1
    fi

    # ========================================================
    # COLLECT DAS TOOL BINS
    # ========================================================

    local BINS=("$BIN_DIR"/*.fa)

    if (( ${#BINS[@]} == 0 )); then

        log "ERROR: DAS Tool bin directory exists but contains no .fa bins"
        log "Upstream DAS Tool output may be incomplete"

        return 1
    fi

    log "DAS Tool bins detected: ${#BINS[@]}"

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then

        log "Removing previous or incomplete CheckM output:"
        log "$OUTDIR"

        rm -rf "$OUTDIR"
    fi

    mkdir -p "$INPUT_DIR"

    rm -f "$LOGFILE"

    # ========================================================
    # FILTER SMALL BINS
    # ========================================================

    log "Filtering bins smaller than ${MIN_SIZE} bp"

    for BIN in "${BINS[@]}"; do

        [[ -s "$BIN" ]] || continue

        # Calculate total sequence length robustly
        SIZE=$(awk '
            !/^>/ {
                gsub(/[[:space:]]/, "")
                total += length($0)
            }
            END {
                print total + 0
            }
        ' "$BIN")

        BIN_NAME=$(basename "$BIN")

        if (( SIZE >= MIN_SIZE )); then

            cp "$BIN" "$INPUT_DIR/$BIN_NAME"

            N_KEPT=$((N_KEPT + 1))

            log "KEEP: $BIN_NAME (${SIZE} bp)"

        else

            log "DROP: $BIN_NAME (${SIZE} bp < ${MIN_SIZE} bp)"

        fi

    done

    log "Bins retained after size filtering: $N_KEPT"

    # ========================================================
    # NO BINS AFTER SIZE FILTER
    # ========================================================

    if (( N_KEPT == 0 )); then

        log "No DAS Tool bins >= ${MIN_SIZE} bp"
        log "CheckM not required"

        touch "$OUTDIR/NO_BINS"

        return 0
    fi

    # ========================================================
    # RUN CHECKM LINEAGE WORKFLOW
    # ========================================================

    log "Running CheckM lineage_wf"
    log "Input bins: $INPUT_DIR"

    set +e

    checkm lineage_wf \
        -x fa \
        -t "$THREADS" \
        "$INPUT_DIR" \
        "$OUTDIR" \
        &>> "$LOGFILE"

    EXIT_CODE=$?

    set -e

    log "CheckM lineage_wf exit code: $EXIT_CODE"

    if (( EXIT_CODE != 0 )); then

        log "ERROR: CheckM lineage_wf failed for $SAMPLE"
        log "Last 50 lines of CheckM log:"

        tail -n 50 "$LOGFILE" || true

        return "$EXIT_CODE"
    fi

    # ========================================================
    # VALIDATE LINEAGE WORKFLOW
    # ========================================================

    if [[ ! -s "$LINEAGE_FILE" ]]; then
        log "ERROR: CheckM lineage file missing:"
        log "$LINEAGE_FILE"
        return 1
    fi

    if [[ ! -s "$BIN_STATS" ]]; then
        log "ERROR: CheckM bin statistics missing:"
        log "$BIN_STATS"
        return 1
    fi

    # ========================================================
    # GENERATE PERSISTENT QA TABLE
    # ========================================================

    log "Generating CheckM QA table"

    set +e

    checkm qa \
        "$LINEAGE_FILE" \
        "$OUTDIR" \
        -o 2 \
        --tab_table \
        -f "$QA_FILE" \
        &>> "$LOGFILE"

    QA_EXIT=$?

    set -e

    log "CheckM QA exit code: $QA_EXIT"

    if (( QA_EXIT != 0 )); then

        log "ERROR: CheckM QA failed for $SAMPLE"
        log "Last 50 lines of CheckM log:"

        tail -n 50 "$LOGFILE" || true

        return "$QA_EXIT"
    fi

    # ========================================================
    # FINAL VALIDATION
    # ========================================================

    if [[ ! -s "$QA_FILE" ]]; then

        log "ERROR: CheckM QA table was not created:"
        log "$QA_FILE"

        return 1
    fi

    if [[ ! -s "$LINEAGE_FILE" ]]; then

        log "ERROR: CheckM lineage file missing"

        return 1
    fi

    if [[ ! -s "$BIN_STATS" ]]; then

        log "ERROR: CheckM bin statistics missing"

        return 1
    fi

    # ========================================================
    # COMPLETE MARKER
    # ========================================================

    touch "$OUTDIR/CHECKM_COMPLETE"

    log "SUCCESS: CheckM completed for $SAMPLE"
    log "QA table:"
    log "$QA_FILE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"
    log "Sample: $SAMPLE"

    run_checkm "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    #
    # IMPORTANT:
    # Use a stable upstream source for array indexing.
    # Do NOT generate the array from only samples currently
    # containing DAS Tool bins.
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$ASM_CLEAN_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.clean.fa" \
            -size +0c \
            | sort
    )

    N=${#ASM_FILES[@]}

    log "Total cleaned assemblies detected: $N"

    if (( N == 0 )); then
        log "ERROR: No cleaned assemblies found"
        exit 1
    fi

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$ASM" .contigs.clean.fa)

    log "Array task ${SLURM_ARRAY_TASK_ID} / $N -> $SAMPLE"

    run_checkm "$SAMPLE"

fi

log "CheckM prokaryotic quality assessment step completed"

;;


# ============================================================
### CHECKM CLASSIFICATION
### FROM DAS TOOL BINS — GLOBAL / NON-ARRAY STEP
# ============================================================
checkm_classify_prok)

set -euo pipefail
shopt -s nullglob

log "Starting CheckM-based MAG classification"

# ============================================================
# CONFIG
# ============================================================

CHECKM_DIR="$OUTPUT_DIR/checkm_prok"
DAS_DIR="$OUTPUT_DIR/das_tool_prok"

CLASS_DIR="$OUTPUT_DIR/checkm_classified_prok"

# Temporary working directory.
# Existing successful classification is preserved until
# the new classification finishes successfully.
WORK_DIR="${CLASS_DIR}.tmp.$$"

SUMMARY="$WORK_DIR/checkm_summary.tsv"

# ============================================================
# INPUT CHECKS
# ============================================================

[[ -d "$CHECKM_DIR" ]] || {
    log "ERROR: CheckM directory not found:"
    log "$CHECKM_DIR"
    exit 1
}

[[ -d "$DAS_DIR" ]] || {
    log "ERROR: DAS Tool directory not found:"
    log "$DAS_DIR"
    exit 1
}

# ============================================================
# PREP TEMPORARY OUTPUT
# ============================================================

rm -rf "$WORK_DIR"

mkdir -p \
    "$WORK_DIR/HQ" \
    "$WORK_DIR/MQ" \
    "$WORK_DIR/LQ"

# ============================================================
# SUMMARY HEADER
# ============================================================

printf "Sample\tBin\tCompleteness\tContamination\tCategory\n" \
    > "$SUMMARY"

# ============================================================
# COUNTERS
# ============================================================

N_SAMPLES=0
N_BINS_TOTAL=0
N_MISSING_BIN=0

# ============================================================
# LOOP OVER CHECKM SAMPLE DIRECTORIES
# ============================================================

for SAMPLE_DIR in "$CHECKM_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")

    QA_FILE="$SAMPLE_DIR/${SAMPLE}_qa.tsv"

    DAS_SAMPLE_DIR="$DAS_DIR/$SAMPLE"
    BIN_DIR="$DAS_SAMPLE_DIR/${SAMPLE}_DASTool_bins"

    # --------------------------------------------------------
    # VALID NO-BIN SAMPLE
    # --------------------------------------------------------

    if [[ -f "$SAMPLE_DIR/NO_BINS" ]]; then
        log "$SAMPLE: CheckM reports NO_BINS -> skipping"
        continue
    fi

    # --------------------------------------------------------
    # QA INPUT
    # --------------------------------------------------------

    if [[ ! -s "$QA_FILE" ]]; then
        log "WARNING: Missing CheckM QA file for $SAMPLE -> skipping"
        continue
    fi

    if [[ ! -d "$BIN_DIR" ]]; then
        log "WARNING: Missing DAS Tool bin directory for $SAMPLE -> skipping"
        continue
    fi

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "QA file: $QA_FILE"
    log "Bin directory: $BIN_DIR"

    TMP_MAP="$WORK_DIR/${SAMPLE}.classification.tmp.tsv"

    rm -f "$TMP_MAP"

    # ========================================================
    # PARSE CHECKM QA TABLE
    #
    # QA was generated previously with:
    #
    #   checkm qa ... --tab_table -o 2
    #
    # Therefore parse by column NAME instead of assuming
    # fixed numeric positions.
    # ========================================================

    awk \
        -F '\t' \
        -v OFS='\t' \
        -v sample="$SAMPLE" \
        -v summary="$SUMMARY" \
        -v map="$TMP_MAP" \
    '
    BEGIN {
        comp_col = 0
        cont_col = 0
        bin_col  = 1
    }

    # --------------------------------------------------------
    # Detect header and required columns
    # --------------------------------------------------------
    /Bin Id/ {

        for (i = 1; i <= NF; i++) {

            field = $i

            gsub(/^[[:space:]]+/, "", field)
            gsub(/[[:space:]]+$/, "", field)

            if (field == "Bin Id")
                bin_col = i

            if (field == "Completeness")
                comp_col = i

            if (field == "Contamination")
                cont_col = i
        }

        next
    }

    # --------------------------------------------------------
    # Skip non-data lines
    # --------------------------------------------------------
    /^[[:space:]]*$/ {
        next
    }

    /^\[/ {
        next
    }

    /^-+/ {
        next
    }

    # --------------------------------------------------------
    # Data rows
    # --------------------------------------------------------
    {

        if (comp_col == 0 || cont_col == 0)
            next

        bin  = $bin_col
        comp = $comp_col
        cont = $cont_col

        # Trim whitespace
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", bin)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", comp)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", cont)

        # Require numeric completeness and contamination
        if (comp !~ /^[0-9]+([.][0-9]+)?$/)
            next

        if (cont !~ /^[0-9]+([.][0-9]+)?$/)
            next

        comp += 0
        cont += 0

        # ----------------------------------------------------
        # CheckM-based quality classification
        # ----------------------------------------------------

        if (comp >= 90 && cont <= 5)
            category = "HQ"

        else if (comp >= 50 && cont <= 10)
            category = "MQ"

        else
            category = "LQ"

        # ----------------------------------------------------
        # Global summary
        # ----------------------------------------------------

        print sample, bin, comp, cont, category >> summary

        # ----------------------------------------------------
        # Bin -> category mapping
        # ----------------------------------------------------

        print bin, category >> map
    }

    END {

        if (comp_col == 0) {
            print "ERROR: Completeness column not found" > "/dev/stderr"
            exit 2
        }

        if (cont_col == 0) {
            print "ERROR: Contamination column not found" > "/dev/stderr"
            exit 2
        }
    }
    ' "$QA_FILE"

    # ========================================================
    # MAKE SURE QA CONTAINED BIN RESULTS
    # ========================================================

    if [[ ! -s "$TMP_MAP" ]]; then

        log "ERROR: No CheckM bin records parsed for $SAMPLE"
        log "Check QA file:"
        log "$QA_FILE"

        exit 1
    fi

    N_SAMPLES=$((N_SAMPLES + 1))

    # ========================================================
    # COPY CLASSIFIED BINS
    # ========================================================

    while IFS=$'\t' read -r BIN CATEGORY; do

        [[ -n "$BIN" ]] || continue
        [[ -n "$CATEGORY" ]] || continue

        SOURCE_BIN="$BIN_DIR/${BIN}.fa"

        if [[ -s "$SOURCE_BIN" ]]; then

            DEST_BIN="$WORK_DIR/$CATEGORY/${SAMPLE}_${BIN}.fa"

            cp "$SOURCE_BIN" "$DEST_BIN"

            N_BINS_TOTAL=$((N_BINS_TOTAL + 1))

        else

            log "WARNING: CheckM bin not found in DAS Tool output:"
            log "$SOURCE_BIN"

            N_MISSING_BIN=$((N_MISSING_BIN + 1))
        fi

    done < "$TMP_MAP"

    rm -f "$TMP_MAP"

done

# ============================================================
# VERIFY THAT SOMETHING WAS CLASSIFIED
# ============================================================

if (( N_SAMPLES == 0 )); then

    log "ERROR: No samples with usable CheckM QA results were found"

    rm -rf "$WORK_DIR"

    exit 1
fi

# ============================================================
# FINAL COUNTS
# ============================================================

HQ_N=$(find "$WORK_DIR/HQ" \
    -maxdepth 1 \
    -type f \
    -name "*.fa" \
    | wc -l)

MQ_N=$(find "$WORK_DIR/MQ" \
    -maxdepth 1 \
    -type f \
    -name "*.fa" \
    | wc -l)

LQ_N=$(find "$WORK_DIR/LQ" \
    -maxdepth 1 \
    -type f \
    -name "*.fa" \
    | wc -l)

SUMMARY_N=$(awk 'NR > 1 {n++} END {print n+0}' "$SUMMARY")

log "----------------------------------------"
log "Samples classified: $N_SAMPLES"
log "CheckM records: $SUMMARY_N"
log "Bins copied: $N_BINS_TOTAL"
log "HQ bins: $HQ_N"
log "MQ bins: $MQ_N"
log "LQ bins: $LQ_N"
log "Missing source bins: $N_MISSING_BIN"

# ============================================================
# CONSISTENCY CHECK
# ============================================================

COPIED_N=$((HQ_N + MQ_N + LQ_N))

if (( COPIED_N != N_BINS_TOTAL )); then

    log "ERROR: Internal bin count mismatch"
    log "Copied counter: $N_BINS_TOTAL"
    log "Files detected: $COPIED_N"

    rm -rf "$WORK_DIR"

    exit 1
fi

if (( N_MISSING_BIN > 0 )); then

    log "ERROR: $N_MISSING_BIN CheckM bins could not be matched"
    log "to the original DAS Tool FASTA files."
    log "Classification will NOT replace the existing output."

    rm -rf "$WORK_DIR"

    exit 1
fi

# ============================================================
# MARK COMPLETE
# ============================================================

touch "$WORK_DIR/CLASSIFICATION_COMPLETE"

# ============================================================
# REPLACE OLD CLASSIFICATION ONLY AFTER SUCCESS
# ============================================================

log "Classification validated successfully"

if [[ -d "$CLASS_DIR" ]]; then

    log "Removing previous classification directory:"
    log "$CLASS_DIR"

    rm -rf "$CLASS_DIR"
fi

mv "$WORK_DIR" "$CLASS_DIR"

# ============================================================
# FINAL OUTPUT
# ============================================================

log "----------------------------------------"
log "CheckM classification completed successfully"
log "Summary:"
log "$CLASS_DIR/checkm_summary.tsv"
log "HQ bins: $HQ_N"
log "MQ bins: $MQ_N"
log "LQ bins: $LQ_N"

;;

# ============================================================
### CHECKM SUMMARY PER SAMPLE
### GLOBAL / NON-ARRAY STEP
# ============================================================
checkm_summary_sample_prok)

set -euo pipefail
shopt -s nullglob

log "Creating CheckM per-sample summary"

# ============================================================
# PATHS
# ============================================================

CLASS_DIR="$OUTPUT_DIR/checkm_classified_prok"

IN="$CLASS_DIR/checkm_summary.tsv"

CHECKM_DIR="$OUTPUT_DIR/checkm_prok"
ASM_CLEAN_DIR="$OUTPUT_DIR/assembly_clean"

OUT="$CLASS_DIR/checkm_summary_per_sample.tsv"

TMP_SAMPLES="$CLASS_DIR/.checkm_samples.$$"
TMP_STATS="$CLASS_DIR/.checkm_stats.$$"

# ============================================================
# INPUT CHECKS
# ============================================================

[[ -s "$IN" ]] || {
    log "ERROR: CheckM classification summary missing:"
    log "$IN"
    exit 1
}

[[ -d "$CHECKM_DIR" ]] || {
    log "ERROR: CheckM directory missing:"
    log "$CHECKM_DIR"
    exit 1
}

[[ -d "$ASM_CLEAN_DIR" ]] || {
    log "ERROR: Clean assembly directory missing:"
    log "$ASM_CLEAN_DIR"
    exit 1
}

# ============================================================
# CLEANUP TEMP FILES ON EXIT
# ============================================================

cleanup_checkm_summary() {
    rm -f "$TMP_SAMPLES" "$TMP_STATS"
}

trap cleanup_checkm_summary EXIT

# ============================================================
# GET STABLE LIST OF ALL SAMPLES
#
# Use assembly_clean rather than checkm_prok so samples that
# produced no CheckM results are still represented.
# ============================================================

find "$ASM_CLEAN_DIR" \
    -mindepth 1 \
    -maxdepth 1 \
    -type d \
    -printf '%f\n' \
    | sort -V \
    > "$TMP_SAMPLES"

N_SAMPLES=$(wc -l < "$TMP_SAMPLES")

if (( N_SAMPLES == 0 )); then
    log "ERROR: No samples found in:"
    log "$ASM_CLEAN_DIR"
    exit 1
fi

log "Total samples detected: $N_SAMPLES"

# ============================================================
# AGGREGATE CHECKM RESULTS
# ============================================================

awk '
BEGIN {
    FS = OFS = "\t"
}

NR == 1 {
    next
}

{
    sample = $1
    comp   = $3 + 0
    cont   = $4 + 0
    cat    = $5

    sum_comp[sample] += comp
    sum_cont[sample] += cont

    n[sample]++

    if (cat == "HQ")
        HQ[sample]++

    else if (cat == "MQ")
        MQ[sample]++

    else if (cat == "LQ")
        LQ[sample]++
}

END {

    for (sample in n) {

        mean_comp = sum_comp[sample] / n[sample]
        mean_cont = sum_cont[sample] / n[sample]

        print sample,
              mean_comp,
              mean_cont,
              n[sample],
              (sample in HQ ? HQ[sample] : 0),
              (sample in MQ ? MQ[sample] : 0),
              (sample in LQ ? LQ[sample] : 0)
    }
}
' "$IN" > "$TMP_STATS"

# ============================================================
# BUILD FINAL TABLE
# ============================================================

printf \
"sample\tmean_completeness\tmean_contamination\tn_bins_assessed\tHQ\tMQ\tLQ\tstatus\n" \
> "$OUT"

awk \
    -F '\t' \
    -v OFS='\t' \
    -v checkm_dir="$CHECKM_DIR" \
'
# ------------------------------------------------------------
# First file = aggregated statistics
# ------------------------------------------------------------

NR == FNR {

    stats[$1] = $0

    next
}

# ------------------------------------------------------------
# Second file = complete stable sample list
# ------------------------------------------------------------

{
    sample = $1

    if (sample in stats) {

        print stats[sample], "classified"

    } else {

        no_bins = checkm_dir "/" sample "/NO_BINS"

        cmd = "test -f \"" no_bins "\""

        if (system(cmd) == 0)
            status = "NO_BINS"
        else
            status = "NO_CHECKM_RESULT"

        print sample,
              "NA",
              "NA",
              0,
              0,
              0,
              0,
              status
    }
}
' "$TMP_STATS" "$TMP_SAMPLES" \
>> "$OUT"

# ============================================================
# VALIDATION
# ============================================================

[[ -s "$OUT" ]] || {
    log "ERROR: Per-sample summary was not created"
    exit 1
}

N_ROWS=$(awk 'NR > 1 {n++} END {print n+0}' "$OUT")

if (( N_ROWS != N_SAMPLES )); then

    log "ERROR: Output row count does not match sample count"
    log "Expected: $N_SAMPLES"
    log "Observed: $N_ROWS"

    exit 1
fi

# ============================================================
# REPORT
# ============================================================

N_CLASSIFIED=$(awk -F'\t' '
    NR > 1 && $8 == "classified" {n++}
    END {print n+0}
' "$OUT")

N_NO_BINS=$(awk -F'\t' '
    NR > 1 && $8 == "NO_BINS" {n++}
    END {print n+0}
' "$OUT")

N_NO_RESULT=$(awk -F'\t' '
    NR > 1 && $8 == "NO_CHECKM_RESULT" {n++}
    END {print n+0}
' "$OUT")

TOTAL_HQ=$(awk -F'\t' '
    NR > 1 {n += $5}
    END {print n+0}
' "$OUT")

TOTAL_MQ=$(awk -F'\t' '
    NR > 1 {n += $6}
    END {print n+0}
' "$OUT")

TOTAL_LQ=$(awk -F'\t' '
    NR > 1 {n += $7}
    END {print n+0}
' "$OUT")

log "----------------------------------------"
log "Per-sample CheckM summary created:"
log "$OUT"
log "Samples represented: $N_ROWS"
log "Samples with classified bins: $N_CLASSIFIED"
log "Samples with NO_BINS: $N_NO_BINS"
log "Samples with no CheckM result: $N_NO_RESULT"
log "Total HQ bins: $TOTAL_HQ"
log "Total MQ bins: $TOTAL_MQ"
log "Total LQ bins: $TOTAL_LQ"

;;

# ============================================================
### BAT (BIN-LEVEL — PROKARYOTIC MAGs)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
bat_prok)

set -euo pipefail
shopt -s nullglob

log "Starting BAT (prokaryotic MAG taxonomy)"

# ============================================================
# CONFIG
# ============================================================

THREADS="${THREADS:-8}"

CLASS_DIR="$OUTPUT_DIR/checkm_classified_prok"
SUMMARY="$CLASS_DIR/checkm_summary.tsv"

BAT_DIR="$OUTPUT_DIR/bat_prok"
LOG_DIR="$OUTPUT_DIR/logs/bat_prok"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$BAT_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "CLASS_DIR: $CLASS_DIR"
log "BAT_DIR: $BAT_DIR"

# ============================================================
# INPUT CHECKS
# ============================================================

[[ -d "$CLASS_DIR" ]] || {
    log "ERROR: CheckM classification directory not found:"
    log "$CLASS_DIR"
    exit 1
}

[[ -s "$SUMMARY" ]] || {
    log "ERROR: CheckM classification summary not found:"
    log "$SUMMARY"
    exit 1
}

[[ -f "$CLASS_DIR/CLASSIFICATION_COMPLETE" ]] || {
    log "ERROR: CheckM classification is not marked complete:"
    log "$CLASS_DIR/CLASSIFICATION_COMPLETE"
    exit 1
}

for CATEGORY in HQ MQ LQ; do

    [[ -d "$CLASS_DIR/$CATEGORY" ]] || {
        log "ERROR: Missing classification directory:"
        log "$CLASS_DIR/$CATEGORY"
        exit 1
    }

done

# ============================================================
# DEPENDENCY CHECKS
# ============================================================

command -v "$CAT_CMD" >/dev/null 2>&1 || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    log "ERROR: DIAMOND not found"
    exit 1
}

command -v python >/dev/null 2>&1 || {
    log "ERROR: Python not found"
    exit 1
}

python -c "import pyrodigal" >/dev/null 2>&1 || {
    log "ERROR: Pyrodigal not available"
    exit 1
}

# ============================================================
# DATABASE CHECKS
# ============================================================

[[ -d "$CAT_DB" ]] || {
    log "ERROR: CAT database not found:"
    log "$CAT_DB"
    exit 1
}

[[ -d "$CAT_TAX" ]] || {
    log "ERROR: CAT taxonomy directory not found:"
    log "$CAT_TAX"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_bat_prok() {

    local SAMPLE="$1"

    local OUTDIR="$BAT_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}"
    local BINS_TMP="$OUTDIR/bins"

    local CLASS_FILE="${PREFIX}.bin2classification.txt"
    local ORF_FILE="${PREFIX}.ORF2LCA.txt"

    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local N_BINS=0
    local EXIT_CODE=0

    # ========================================================
    # SAMPLE TEMP DIRECTORY
    # ========================================================

    local TMP_ROOT="${TMPDIR:-/tmp}"
    local SAMPLE_TMP="${TMP_ROOT}/bat_${SAMPLE}"

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"

    # ========================================================
    # SKIP COMPLETED
    # ========================================================

    if [[ "${FORCE:-0}" != "1" ]]; then

        if [[ -f "$OUTDIR/BAT_COMPLETE" ]]; then
            log "BAT already completed -> skipping"
            return 0
        fi

        if [[ -f "$OUTDIR/NO_BINS" ]]; then
            log "No classified MAGs for this sample -> skipping"
            return 0
        fi

        # Backward compatibility with older BAT results
        if [[ -s "$CLASS_FILE" && -s "$ORF_FILE" ]]; then

            log "Existing complete BAT outputs detected -> skipping"

            touch "$OUTDIR/BAT_COMPLETE"

            return 0
        fi

    fi

    # ========================================================
    # FIND CLASSIFIED MAGS
    # ========================================================

    local SAMPLE_BINS=()

    mapfile -t SAMPLE_BINS < <(
        find \
            "$CLASS_DIR/HQ" \
            "$CLASS_DIR/MQ" \
            "$CLASS_DIR/LQ" \
            -maxdepth 1 \
            -type f \
            -name "${SAMPLE}_*.fa" \
            -size +0c \
            | sort -V
    )

    N_BINS=${#SAMPLE_BINS[@]}

    log "Classified MAGs detected: $N_BINS"

    # ========================================================
    # VALID ZERO-BIN RESULT
    # ========================================================

    if (( N_BINS == 0 )); then

        log "No classified MAGs found for $SAMPLE"
        log "BAT not required"

        rm -rf "$OUTDIR"
        mkdir -p "$OUTDIR"

        touch "$OUTDIR/NO_BINS"

        return 0
    fi

    # ========================================================
    # CLEAN PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then

        log "Removing previous or incomplete BAT output:"
        log "$OUTDIR"

        rm -rf "$OUTDIR"
    fi

    mkdir -p "$BINS_TMP"

    rm -f "$LOGFILE"

    # ========================================================
    # PREPARE SAMPLE MAG DIRECTORY
    #
    # Use symlinks rather than duplicating all MAG FASTAs.
    # ========================================================

    log "Preparing BAT input directory"

    for BIN in "${SAMPLE_BINS[@]}"; do

        [[ -s "$BIN" ]] || continue

        BASENAME=$(basename "$BIN")

        if [[ -e "$BINS_TMP/$BASENAME" ]]; then

            log "ERROR: Duplicate MAG filename detected:"
            log "$BASENAME"

            return 1
        fi

        ln -s "$BIN" "$BINS_TMP/$BASENAME"

    done

    # ========================================================
    # VERIFY BIN COUNT
    # ========================================================

    N_LINKED=$(find "$BINS_TMP" \
        -maxdepth 1 \
        -type l \
        -name "*.fa" \
        | wc -l)

    log "MAGs prepared for BAT: $N_LINKED"

    if (( N_LINKED != N_BINS )); then

        log "ERROR: Input MAG count mismatch"
        log "Expected: $N_BINS"
        log "Prepared: $N_LINKED"

        return 1
    fi

    # ========================================================
    # PREP TEMP DIRECTORY
    # ========================================================

    rm -rf "$SAMPLE_TMP"
    mkdir -p "$SAMPLE_TMP"

    export TMPDIR="$SAMPLE_TMP"
    export TEMP="$SAMPLE_TMP"
    export TMP="$SAMPLE_TMP"

    log "BAT TMPDIR: $TMPDIR"

    # ========================================================
    # RUN BAT
    # ========================================================

    log "Running CAT_pack bins (BAT)"
    log "Input MAGs: $N_BINS"

    set +e

    "$CAT_CMD" bins \
        -b "$BINS_TMP" \
        -d "$CAT_DB" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --bin_suffix .fa \
        --nproc "$THREADS" \
        --block_size 2 \
        --tmpdir "$TMPDIR" \
        2>&1 | tee "$LOGFILE"

    EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "BAT exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if (( EXIT_CODE != 0 )); then

        log "ERROR: BAT failed for $SAMPLE"
        log "Last 50 lines of BAT log:"

        tail -n 50 "$LOGFILE" || true

        return "$EXIT_CODE"
    fi

    # ========================================================
    # VALIDATE OUTPUTS
    # ========================================================

    if [[ ! -s "$CLASS_FILE" ]]; then

        log "ERROR: BAT bin2classification output missing:"
        log "$CLASS_FILE"

        return 1
    fi

    if [[ ! -s "$ORF_FILE" ]]; then

        log "ERROR: BAT ORF2LCA output missing:"
        log "$ORF_FILE"

        return 1
    fi

    # ========================================================
    # MARK COMPLETE
    # ========================================================

    rm -f "$OUTDIR/NO_BINS"

    touch "$OUTDIR/BAT_COMPLETE"

    # ========================================================
    # CLEAN TEMP
    # ========================================================

    rm -rf "$SAMPLE_TMP"

    log "SUCCESS: BAT completed for $SAMPLE"
    log "Classification:"
    log "$CLASS_FILE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"
    log "Sample: $SAMPLE"

    run_bat_prok "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    #
    # Use the finished CheckM classification summary as the
    # stable source of samples rather than discovering samples
    # from changing output directories.
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t SAMPLES < <(
        awk -F '\t' '
            NR > 1 && $1 != "" {
                print $1
            }
        ' "$SUMMARY" \
        | sort -Vu
    )

    N=${#SAMPLES[@]}

    log "Samples with classified MAGs: $N"

    if (( N == 0 )); then
        log "ERROR: No classified MAG samples found"
        exit 1
    fi

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 ||
          SLURM_ARRAY_TASK_ID > N )); then

        log "ERROR: Array task ${SLURM_ARRAY_TASK_ID} outside range 1-$N"

        exit 1
    fi

    SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

    log "Array task ${SLURM_ARRAY_TASK_ID} / $N -> $SAMPLE"

    run_bat_prok "$SAMPLE"

fi

log "BAT prokaryotic MAG taxonomy step completed"

;;

# ============================================================
### GTDB-TK (PROK TAXONOMY — PER SAMPLE)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
gtdbtk)

set -euo pipefail
shopt -s nullglob

log "Starting GTDB-Tk classification"

# ============================================================
# CONFIG
# ============================================================

THREADS="${THREADS:-8}"
PPLACER_CPUS="${PPLACER_CPUS:-4}"

# Never allow pplacer to request more CPUs than the job
if (( PPLACER_CPUS > THREADS )); then
    PPLACER_CPUS="$THREADS"
fi

export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CLASS_DIR="$OUTPUT_DIR/checkm_classified_prok"
ASM_CLEAN_DIR="$OUTPUT_DIR/assembly_clean"

GTDB_DIR="$OUTPUT_DIR/gtdbtk"
LOG_DIR="$OUTPUT_DIR/logs/gtdbtk"

# ============================================================
# QUALITY FILTER
#
# Recommended:
#   HQ + MQ = primary MAG taxonomy
#
# Set INCLUDE_LQ=1 if you intentionally want GTDB-Tk
# classifications for LQ bins as well.
# ============================================================

INCLUDE_LQ="${INCLUDE_LQ:-0}"

# ============================================================
# GTDB DATABASE
# ============================================================

export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-/staging/leuven/stg_00100/Rodgee/databases/gtdbtk/release232}"

mkdir -p "$GTDB_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "PPLACER_CPUS: $PPLACER_CPUS"
log "GTDBTK_DATA_PATH: $GTDBTK_DATA_PATH"
log "INCLUDE_LQ: $INCLUDE_LQ"

# ============================================================
# INPUT CHECKS
# ============================================================

[[ -d "$CLASS_DIR" ]] || {
    log "ERROR: CheckM classification directory not found:"
    log "$CLASS_DIR"
    exit 1
}

[[ -f "$CLASS_DIR/CLASSIFICATION_COMPLETE" ]] || {
    log "ERROR: CheckM classification is not marked complete:"
    log "$CLASS_DIR/CLASSIFICATION_COMPLETE"
    exit 1
}

[[ -d "$ASM_CLEAN_DIR" ]] || {
    log "ERROR: Clean assembly directory not found:"
    log "$ASM_CLEAN_DIR"
    exit 1
}

[[ -d "$GTDBTK_DATA_PATH" ]] || {
    log "ERROR: GTDBTK_DATA_PATH not found:"
    log "$GTDBTK_DATA_PATH"
    exit 1
}

# ============================================================
# TOOL CHECK
# ============================================================

command -v gtdbtk >/dev/null 2>&1 || {
    log "ERROR: gtdbtk not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_gtdbtk() {

    local SAMPLE="$1"

    local OUTDIR="$GTDB_DIR/$SAMPLE"
    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    # GTDB-Tk classify_wf puts final summary files here
    local BAC="$OUTDIR/classify/gtdbtk.bac120.summary.tsv"
    local ARC="$OUTDIR/classify/gtdbtk.ar53.summary.tsv"

    local N_MAGS=0
    local N_REPORTED=0
    local EXIT_CODE=0

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"

    # ========================================================
    # SKIP COMPLETED
    # ========================================================

    if [[ "${FORCE:-0}" != "1" ]]; then

        if [[ -f "$OUTDIR/GTDBTK_COMPLETE" ]]; then
            log "GTDB-Tk already completed -> skipping"
            return 0
        fi

        if [[ -f "$OUTDIR/NO_BINS" ]]; then
            log "No MAGs available for GTDB-Tk -> skipping"
            return 0
        fi

        if [[ -f "$OUTDIR/NO_ELIGIBLE_MAGS" ]]; then
            log "No HQ/MQ MAGs available for GTDB-Tk -> skipping"
            return 0
        fi

        # Backward compatibility
        if [[ -s "$BAC" || -s "$ARC" ]]; then

            log "Existing GTDB-Tk summary detected -> marking complete"

            touch "$OUTDIR/GTDBTK_COMPLETE"

            return 0
        fi

    fi

    # ========================================================
    # DETERMINE QUALITY DIRECTORIES
    # ========================================================

    local QUALITY_DIRS=(
        "$CLASS_DIR/HQ"
        "$CLASS_DIR/MQ"
    )

    if [[ "$INCLUDE_LQ" == "1" ]]; then
        QUALITY_DIRS+=("$CLASS_DIR/LQ")
    fi

    # ========================================================
    # COLLECT SAMPLE MAGs
    # ========================================================

    local SAMPLE_MAGS=()

    mapfile -t SAMPLE_MAGS < <(
        find "${QUALITY_DIRS[@]}" \
            -maxdepth 1 \
            -type f \
            -name "${SAMPLE}_*.fa" \
            -size +0c \
            2>/dev/null \
            | sort -V
    )

    N_MAGS=${#SAMPLE_MAGS[@]}

    log "Eligible MAGs detected: $N_MAGS"

    # ========================================================
    # DISTINGUISH TRUE NO-BIN FROM QUALITY FILTER
    # ========================================================

    if (( N_MAGS == 0 )); then

        local N_ALL=0

        N_ALL=$(find \
            "$CLASS_DIR/HQ" \
            "$CLASS_DIR/MQ" \
            "$CLASS_DIR/LQ" \
            -maxdepth 1 \
            -type f \
            -name "${SAMPLE}_*.fa" \
            -size +0c \
            2>/dev/null \
            | wc -l)

        rm -rf "$OUTDIR"
        mkdir -p "$OUTDIR"

        if (( N_ALL == 0 )); then

            log "No classified bins exist for $SAMPLE"
            log "GTDB-Tk not required"

            touch "$OUTDIR/NO_BINS"

        else

            log "Bins exist, but none pass the selected GTDB-Tk quality filter"
            log "Total classified bins: $N_ALL"

            touch "$OUTDIR/NO_ELIGIBLE_MAGS"

        fi

        return 0
    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then

        log "Removing previous or incomplete GTDB-Tk output:"
        log "$OUTDIR"

        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    rm -f "$LOGFILE"

    # ========================================================
    # SAMPLE-SPECIFIC TEMP AREA
    #
    # Keep input FASTAs OUTSIDE the GTDB-Tk output directory.
    # ========================================================

    local TMP_ROOT="${TMPDIR:-/tmp}"

    local SAMPLE_TMP="$TMP_ROOT/gtdbtk_${SAMPLE}"
    local INPUT_DIR="$SAMPLE_TMP/input"
    local GTDB_TMP="$SAMPLE_TMP/tmp"

    rm -rf "$SAMPLE_TMP"

    mkdir -p \
        "$INPUT_DIR" \
        "$GTDB_TMP"

    local QUALITY_MAP="$OUTDIR/${SAMPLE}.gtdbtk_input_quality.tsv"

    printf "genome\tquality\tsource_file\n" > "$QUALITY_MAP"

    # ========================================================
    # PREPARE INPUT MAGs
    #
    # Preserve original MAG identifiers.
    # Do NOT append _HQ/_MQ to the genome identifier.
    # ========================================================

    log "Preparing GTDB-Tk MAG input"

    for MAG in "${SAMPLE_MAGS[@]}"; do

        [[ -s "$MAG" ]] || continue

        local BASENAME
        local GENOME
        local QUALITY

        BASENAME=$(basename "$MAG")
        GENOME="${BASENAME%.fa}"
        QUALITY=$(basename "$(dirname "$MAG")")

        if [[ -e "$INPUT_DIR/$BASENAME" ]]; then

            log "ERROR: Duplicate MAG filename:"
            log "$BASENAME"

            return 1
        fi

        ln -s "$MAG" "$INPUT_DIR/$BASENAME"

        printf "%s\t%s\t%s\n" \
            "$GENOME" \
            "$QUALITY" \
            "$MAG" \
            >> "$QUALITY_MAP"
    done

    # ========================================================
    # VERIFY INPUT COUNT
    # ========================================================

    local N_LINKED=0

    N_LINKED=$(find "$INPUT_DIR" \
        -maxdepth 1 \
        -type l \
        -name "*.fa" \
        | wc -l)

    log "GTDB-Tk input MAGs prepared: $N_LINKED"

    if (( N_LINKED != N_MAGS )); then

        log "ERROR: MAG input count mismatch"
        log "Expected: $N_MAGS"
        log "Prepared: $N_LINKED"

        return 1
    fi

    # ========================================================
    # RUN GTDB-TK
    # ========================================================

    log "Running GTDB-Tk classify_wf"
    log "MAGs: $N_MAGS"
    log "CPUs: $THREADS"
    log "pplacer CPUs: $PPLACER_CPUS"
    log "Temporary directory: $GTDB_TMP"

    set +e

    TMPDIR="$GTDB_TMP" \
    TEMP="$GTDB_TMP" \
    TMP="$GTDB_TMP" \
    gtdbtk classify_wf \
        --genome_dir "$INPUT_DIR" \
        --out_dir "$OUTDIR" \
        --extension fa \
        --cpus "$THREADS" \
        --pplacer_cpus "$PPLACER_CPUS" \
        --tmpdir "$GTDB_TMP" \
        &>> "$LOGFILE"

    EXIT_CODE=$?

    set -e

    log "GTDB-Tk exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if (( EXIT_CODE != 0 )); then

        log "ERROR: GTDB-Tk failed for $SAMPLE"
        log "Last 50 lines of GTDB-Tk log:"

        tail -n 50 "$LOGFILE" || true

        log "Temporary files retained at:"
        log "$SAMPLE_TMP"

        return "$EXIT_CODE"
    fi

    # ========================================================
    # VALIDATE OUTPUT
    # ========================================================

    if [[ ! -s "$BAC" && ! -s "$ARC" ]]; then

        log "ERROR: Neither bacterial nor archaeal summary was produced"
        log "Expected one or both of:"
        log "$BAC"
        log "$ARC"

        return 1
    fi

    # ========================================================
    # COUNT REPORTED GENOMES
    # ========================================================

    N_REPORTED=0

    if [[ -s "$BAC" ]]; then

        N_BAC=$(awk '
            NR > 1 {n++}
            END {print n+0}
        ' "$BAC")

        log "Bacterial genomes reported: $N_BAC"

        N_REPORTED=$((N_REPORTED + N_BAC))
    fi

    if [[ -s "$ARC" ]]; then

        N_ARC=$(awk '
            NR > 1 {n++}
            END {print n+0}
        ' "$ARC")

        log "Archaeal genomes reported: $N_ARC"

        N_REPORTED=$((N_REPORTED + N_ARC))
    fi

    log "Input MAGs: $N_MAGS"
    log "GTDB-Tk summary records: $N_REPORTED"

    # Do not hard-fail solely on count mismatch because
    # GTDB-Tk versions may handle failed/filtered genomes
    # differently. Make it highly visible instead.
    if (( N_REPORTED != N_MAGS )); then

        log "WARNING: GTDB-Tk output count differs from input MAG count"
        log "Check:"
        log "$OUTDIR/identify/gtdbtk.failed_genomes.tsv"
        log "$LOGFILE"
    fi

    # ========================================================
    # MARK COMPLETE
    # ========================================================

    rm -f \
        "$OUTDIR/NO_BINS" \
        "$OUTDIR/NO_ELIGIBLE_MAGS"

    touch "$OUTDIR/GTDBTK_COMPLETE"

    # ========================================================
    # CLEAN TEMP
    # ========================================================

    rm -rf "$SAMPLE_TMP"

    # ========================================================
    # REPORT
    # ========================================================

    if [[ -s "$BAC" ]]; then
        log "Bacterial classifications:"
        log "$BAC"
    fi

    if [[ -s "$ARC" ]]; then
        log "Archaeal classifications:"
        log "$ARC"
    fi

    log "SUCCESS: GTDB-Tk completed for $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"
    log "Sample: $SAMPLE"

    run_gtdbtk "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    #
    # Stable array based on cleaned assemblies.
    # This represents all biological samples, including
    # samples such as S8/S22 that ultimately have no MAGs.
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t ASM_FILES < <(
        find "$ASM_CLEAN_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -name "*.contigs.clean.fa" \
            -size +0c \
            | sort -V
    )

    N=${#ASM_FILES[@]}

    log "Total samples detected: $N"

    if (( N == 0 )); then
        log "ERROR: No cleaned assemblies found"
        exit 1
    fi

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 ||
          SLURM_ARRAY_TASK_ID > N )); then

        log "ERROR: Array task ${SLURM_ARRAY_TASK_ID} outside range 1-$N"

        exit 1
    fi

    ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$ASM" .contigs.clean.fa)

    log "Array task ${SLURM_ARRAY_TASK_ID} / $N -> $SAMPLE"

    run_gtdbtk "$SAMPLE"

fi

log "GTDB-Tk prokaryotic taxonomy step completed"

;;

# ============================================================
### GTDB SUMMARY (MERGED TAXONOMY — ALL MAGS)
# ============================================================
gtdbtk_summary)

set -euo pipefail
shopt -s nullglob

log "Creating GTDB summary"

GTDB_DIR="$OUTPUT_DIR/gtdbtk"

OUT="$GTDB_DIR/gtdbtk_summary.tsv"

echo -e "MAG\tGenome\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > "$OUT"

# ============================================================
# LOOP OVER GTDB SUMMARY FILES
# ============================================================

for FILE in "$GTDB_DIR"/*/classify/*.summary.tsv; do

    [[ -s "$FILE" ]] || continue

    SAMPLE=$(basename "$(dirname "$(dirname "$FILE")")")

    log "Processing: $SAMPLE"

    awk -F'\t' -v sample="$SAMPLE" '
    BEGIN {
        OFS="\t"
    }

    # --------------------------------------------------------
    # SKIP HEADER
    # --------------------------------------------------------
    NR > 1 {

        genome = $1

        taxonomy = $2

        split(taxonomy, tax, ";")

        d="NA"
        p="NA"
        c="NA"
        o="NA"
        f="NA"
        g="NA"
        s="NA"

        for (i in tax) {

            if (tax[i] ~ /^d__/) d=tax[i]

            else if (tax[i] ~ /^p__/) p=tax[i]

            else if (tax[i] ~ /^c__/) c=tax[i]

            else if (tax[i] ~ /^o__/) o=tax[i]

            else if (tax[i] ~ /^f__/) f=tax[i]

            else if (tax[i] ~ /^g__/) g=tax[i]

            else if (tax[i] ~ /^s__/) s=tax[i]
        }

        print sample, genome, d, p, c, o, f, g, s
    }
    ' "$FILE" >> "$OUT"

done

# ============================================================
# REPORT
# ============================================================

N=$(($(wc -l < "$OUT") - 1))

log "GTDB summary created: $OUT"

log "Total MAG classifications: $N"

;;

# ============================================================
### MERGE CHECKM + GTDB (FINAL TABLE)
# ============================================================
gtdbtk_merge_checkm)

set -euo pipefail

CHECKM="$OUTPUT_DIR/checkm_classified_prok/checkm_summary.tsv"

GTDB="$OUTPUT_DIR/gtdbtk/gtdbtk_summary.tsv"

OUT="$OUTPUT_DIR/gtdbtk/gtdbtk_checkm_merged.tsv"

log "Merging CheckM + GTDB"

awk '
BEGIN {
    FS=OFS="\t"
}

# ============================================================
# LOAD GTDB TAXONOMY
# ============================================================
NR==FNR {

    if (FNR==1) next

    genome=$2

    # --------------------------------------------------------
    # CLEAN GTDB GENOME NAME
    # --------------------------------------------------------
    gsub(/\.fa_HQ$/, "", genome)
    gsub(/\.fa_MQ$/, "", genome)
    gsub(/\.fa_LQ$/, "", genome)
    gsub(/\.fa$/, "", genome)

    tax[genome]=$0

    next
}

# ============================================================
# PROCESS CHECKM
# ============================================================
FNR==1 {

    print $0,"Domain","Phylum","Class","Order","Family","Genus","Species"

    next
}

{

    sample=$1
    bin=$2
    # --------------------------------------------------------
    # MATCH THE CLASSIFIED FASTA ID: SAMPLE_BIN
    # --------------------------------------------------------
    key=sample "_" bin

    # Older GTDB-Tk runs appended the CheckM quality to the FASTA ID.
    if (!(key in tax)) key=key "_" $5

    if (key in tax) {

        split(tax[key], t, "\t")

        print $0, t[3], t[4], t[5], t[6], t[7], t[8], t[9]

    } else {

        print $0,"NA","NA","NA","NA","NA","NA","NA"
    }
}
' "$GTDB" "$CHECKM" > "$OUT"

log "Merged file created: $OUT"

;;

# ============================================================
### TIARA SPLIT (PER-SAMPLE; ARRAY OR SINGLE-SAMPLE)
# ============================================================
tiara_split)

set -euo pipefail
shopt -s nullglob

log "Starting Tiara split step"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-8}"

ASM_DIR="$OUTPUT_DIR/assembly"
TIARA_DIR="$OUTPUT_DIR/tiara"

mkdir -p "$TIARA_DIR"

log "THREADS: $THREADS"
log "Assembly directory: $ASM_DIR"
log "Tiara directory: $TIARA_DIR"

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    CONTIGS="$ASM_DIR/$SAMPLE/${SAMPLE}.contigs.fa"

    [[ -s "$CONTIGS" ]] || {
        log "ERROR: Assembly not found:"
        log "$CONTIGS"
        exit 1
    }

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    # -------------------- COLLECT ASSEMBLIES --------------------
    mapfile -t ASM_FILES < <(
        for d in "$ASM_DIR"/*; do

            [[ -d "$d" ]] || continue

            s=$(basename "$d")
            f="$d/${s}.contigs.fa"

            [[ -s "$f" ]] && echo "$f"

        done | sort
    )

    N=${#ASM_FILES[@]}

    log "Total assemblies detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No valid assemblies found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    CONTIGS="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    SAMPLE=$(basename "$CONTIGS" | sed -E 's/\.contigs\.fa$//')

fi

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUT="$TIARA_DIR/$SAMPLE"

mkdir -p "$OUT"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "Contigs: $CONTIGS"
log "Output: $OUT"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# CHECK TIARA
# ============================================================

command -v tiara >/dev/null || {
    log "ERROR: tiara not found"
    log "Make sure the Tiara environment is activated in SLURM"
    exit 1
}

tiara --help >/dev/null 2>&1 || {
    log "ERROR: tiara command found but not working"
    exit 1
}

# ============================================================
# SKIP COMPLETED SAMPLE
# ============================================================

if [[ -s "$OUT/tiara.txt" && "${FORCE:-0}" != "1" ]]; then
    log "Tiara already completed -> skipping"
    log "$OUT/tiara.txt"
    exit 0
fi

# ============================================================
# CLEAN PREVIOUS PARTIAL OUTPUTS
# ============================================================

if [[ "${FORCE:-0}" == "1" ]]; then

    log "FORCE=1 -> cleaning previous Tiara outputs"

    rm -f \
        "$OUT/tiara.txt" \
        "$OUT/euk_ids.txt" \
        "$OUT/prok_ids.txt" \
        "$OUT/organelle_ids.txt" \
        "$OUT/mito_ids.txt" \
        "$OUT/plastid_ids.txt" \
        "$OUT/euk.fa" \
        "$OUT/prok.fa" \
        "$OUT/organelle.fa" \
        "$OUT/mito.fa" \
        "$OUT/plastid.fa"

fi

# ============================================================
# RUN TIARA
# ============================================================

log "Running Tiara"

tiara \
    -i "$CONTIGS" \
    -o "$OUT/tiara.txt" \
    -t "$THREADS" \
    --min_len 3000

# ============================================================
# VALIDATE
# ============================================================

[[ -s "$OUT/tiara.txt" ]] || {
    log "ERROR: Tiara failed"
    exit 1
}

# ============================================================
# SPLIT CONTIG IDS
# ============================================================

log "Extracting Tiara classification groups"

# Eukaryotic
awk 'NR>1 && $5=="eukarya"{print $1}' \
    "$OUT/tiara.txt" > "$OUT/euk_ids.txt"

# Prokaryotic
awk 'NR>1 && ($5=="bacteria" || $5=="archaea" || $5=="prokarya"){print $1}' \
    "$OUT/tiara.txt" > "$OUT/prok_ids.txt"

# All organelles
awk 'NR>1 && $5=="organelle"{print $1}' \
    "$OUT/tiara.txt" > "$OUT/organelle_ids.txt"

# Mitochondria
awk 'NR>1 && $5=="organelle" && $6=="mitochondrion"{print $1}' \
    "$OUT/tiara.txt" > "$OUT/mito_ids.txt"

# Plastids
awk 'NR>1 && $5=="organelle" && $6=="plastid"{print $1}' \
    "$OUT/tiara.txt" > "$OUT/plastid_ids.txt"

# ============================================================
# FASTA SPLIT
# ============================================================

log "Splitting FASTA (single-pass)"

# Remove previous split FASTAs first.
# Important if a failed/partial run is being repeated.
rm -f \
    "$OUT/euk.fa" \
    "$OUT/prok.fa" \
    "$OUT/organelle.fa" \
    "$OUT/mito.fa" \
    "$OUT/plastid.fa"

awk \
-v euk="$OUT/euk_ids.txt" \
-v prok="$OUT/prok_ids.txt" \
-v org="$OUT/organelle_ids.txt" \
-v mito="$OUT/mito_ids.txt" \
-v plast="$OUT/plastid_ids.txt" \
-v euk_out="$OUT/euk.fa" \
-v prok_out="$OUT/prok.fa" \
-v org_out="$OUT/organelle.fa" \
-v mito_out="$OUT/mito.fa" \
-v plast_out="$OUT/plastid.fa" \
'
BEGIN {

    while ((getline < euk) > 0)
        euk_ids[$1]=1

    while ((getline < prok) > 0)
        prok_ids[$1]=1

    while ((getline < org) > 0)
        org_ids[$1]=1

    while ((getline < mito) > 0)
        mito_ids[$1]=1

    while ((getline < plast) > 0)
        plast_ids[$1]=1

    close(euk)
    close(prok)
    close(org)
    close(mito)
    close(plast)
}

/^>/ {

    split(substr($0,2), a, " ")

    header=a[1]

    is_euk   = (header in euk_ids)
    is_prok  = (header in prok_ids)
    is_org   = (header in org_ids)
    is_mito  = (header in mito_ids)
    is_plast = (header in plast_ids)
}

{
    if (is_euk)
        print >> euk_out

    if (is_prok)
        print >> prok_out

    if (is_org)
        print >> org_out

    if (is_mito)
        print >> mito_out

    if (is_plast)
        print >> plast_out
}

END {

    close(euk_out)
    close(prok_out)
    close(org_out)
    close(mito_out)
    close(plast_out)
}
' "$CONTIGS"

# ============================================================
# EMPTY FILE SAFETY
# ============================================================

for f in \
    euk.fa \
    prok.fa \
    organelle.fa \
    mito.fa \
    plastid.fa
do

    [[ -s "$OUT/$f" ]] || : > "$OUT/$f"

done

# ============================================================
# SUMMARY
# ============================================================

log "Summary for $SAMPLE:"

log "  euk:       $(wc -l < "$OUT/euk_ids.txt")"
log "  prok:      $(wc -l < "$OUT/prok_ids.txt")"
log "  organelle: $(wc -l < "$OUT/organelle_ids.txt")"
log "  mito:      $(wc -l < "$OUT/mito_ids.txt")"
log "  plastid:   $(wc -l < "$OUT/plastid_ids.txt")"

log "Finished Tiara split for $SAMPLE"

;;

# ============================================================
### REMAG (EUK BINNING — PER-SAMPLE; ARRAY OR SINGLE-SAMPLE)
# ============================================================
remag)

set -euo pipefail
shopt -s nullglob

log "Starting reMAG step"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-8}"

TIARA_DIR="$OUTPUT_DIR/tiara"
MAP_DIR="$OUTPUT_DIR/mapping"
REMAG_DIR="$OUTPUT_DIR/remag"
LOG_DIR="$OUTPUT_DIR/logs/remag"

mkdir -p "$REMAG_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "TIARA_DIR: $TIARA_DIR"
log "MAP_DIR: $MAP_DIR"
log "REMAG_DIR: $REMAG_DIR"

# ============================================================
# CHECK INPUT DIRECTORIES
# ============================================================

[[ -d "$TIARA_DIR" ]] || {
    log "ERROR: Missing TIARA directory:"
    log "$TIARA_DIR"
    exit 1
}

[[ -d "$MAP_DIR" ]] || {
    log "ERROR: Missing mapping directory:"
    log "$MAP_DIR"
    exit 1
}

# ============================================================
# TOOL CHECK
# ============================================================

command -v remag >/dev/null 2>&1 || {
    log "ERROR: remag not found"
    log "Make sure the reMAG conda environment is activated in SLURM"
    exit 1
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    log "Running in single-sample mode"

    COND="$TIARA_DIR/$SAMPLE"

    [[ -d "$COND" ]] || {
        log "ERROR: Tiara directory not found:"
        log "$COND"
        exit 1
    }

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------
    log "Running in array mode"

    mapfile -t TIARA_SAMPLES < <(
        find "$TIARA_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#TIARA_SAMPLES[@]}

    log "Total Tiara sample directories detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No Tiara sample directories found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------
    COND="${TIARA_SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$COND")

fi

# ============================================================
# SAMPLE FILES
# ============================================================

EUK_FASTA="$COND/euk.fa"
BAM="$MAP_DIR/${SAMPLE}.bam"

OUTDIR="$REMAG_DIR/$SAMPLE"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "EUK_FASTA: $EUK_FASTA"
log "BAM: $BAM"
log "OUTDIR: $OUTDIR"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Array task: ${SLURM_ARRAY_TASK_ID}"
fi

# ============================================================
# INPUT VALIDATION
# ============================================================

if [[ ! -s "$EUK_FASTA" ]]; then
    log "WARNING: Missing or empty euk.fa -> skipping"
    continue
fi

if [[ ! -s "$BAM" ]]; then
    log "WARNING: Missing or empty BAM -> skipping"
    continue
fi

# ============================================================
# SKIP IF DONE
# ============================================================

if [[ -s "$OUTDIR/bins.csv" && "${FORCE:-0}" != "1" ]]; then
    log "reMAG already completed -> skipping"
    log "$OUTDIR/bins.csv"
    continue
fi

# ============================================================
# CLEAN OLD / PARTIAL RUN
# ============================================================

if [[ -d "$OUTDIR" ]]; then
    log "Removing previous reMAG directory:"
    log "$OUTDIR"
    rm -rf "$OUTDIR"
fi

mkdir -p "$OUTDIR"

# ============================================================
# RUN REMAG
# ============================================================

log "Running reMAG"

REMAG_LOG="$LOG_DIR/${SAMPLE}.log"

REMAG_EXIT=0

remag \
    "$EUK_FASTA" \
    -c "$BAM" \
    -o "$OUTDIR" \
    --skip-bacterial-filter \
    &>> "$REMAG_LOG" || REMAG_EXIT=$?

log "reMAG exit code: $REMAG_EXIT"

# ============================================================
# CHECK REMAG RESULT
# ============================================================

# ------------------------------------------------------------
# CASE 1: reMAG produced bins.csv
# ------------------------------------------------------------

if [[ -s "$OUTDIR/bins.csv" ]]; then

    if [[ -d "$OUTDIR/bins" ]]; then

        NBINS=$(find "$OUTDIR/bins" \
            -maxdepth 1 \
            -type f \
            -name "*.fa" \
            | wc -l)

    else
        NBINS=0
    fi

    if (( NBINS > 0 )); then

        log "reMAG completed successfully"
        log "Number of bins: $NBINS"

        touch "$OUTDIR/REMAG_COMPLETE"

    else

        log "reMAG completed but produced no bins"
        log "This is treated as a valid no-bin result"

        touch "$OUTDIR/NO_BINS"

    fi

    continue
fi

# ------------------------------------------------------------
# CASE 2: exit code 0 but no bins.csv
# ------------------------------------------------------------

if (( REMAG_EXIT == 0 )); then

    log "reMAG completed without producing bins.csv"
    log "No bins were recovered for $SAMPLE"

    touch "$OUTDIR/NO_BINS"

    continue
fi

# ------------------------------------------------------------
# CASE 3: genuine reMAG failure
# ------------------------------------------------------------

log "ERROR: reMAG failed for $SAMPLE"
log "Exit code: $REMAG_EXIT"
log "Check log:"
log "$REMAG_LOG"

log "Last 50 lines of reMAG log:"
tail -n 50 "$REMAG_LOG" || true

exit "$REMAG_EXIT"

# ============================================================
# VALIDATE
# ============================================================

if [[ ! -s "$OUTDIR/bins.csv" ]]; then
    log "WARNING: reMAG finished but no bins.csv was produced"
    log "Sample: $SAMPLE"
    log "Check log:"
    log "$LOG_DIR/${SAMPLE}.log"
    continue
fi

# ============================================================
# SUMMARY
# ============================================================

log "Finished reMAG for $SAMPLE"

if [[ -d "$OUTDIR/bins" ]]; then

    NBINS=$(find "$OUTDIR/bins" \
        -maxdepth 1 \
        -type f \
        -name "*.fa" \
        | wc -l)

    log "Number of bins: $NBINS"

else

    log "WARNING: No bins directory created"

fi

log "reMAG completed for sample: $SAMPLE"

;;

# ============================================================
### BUSCO (EUK QUALITY — PER-BIN; ARRAY OR SINGLE-SAMPLE)
# ============================================================
busco)

set -euo pipefail
shopt -s nullglob

log "Starting BUSCO step"

THREADS="${THREADS:-8}"
MIN_SIZE="${BUSCO_MIN_SIZE:-50000}"

REMAG_DIR="$OUTPUT_DIR/remag"
BUSCO_DIR="$OUTPUT_DIR/busco"
LOG_DIR="$OUTPUT_DIR/logs/busco"

mkdir -p "$BUSCO_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "REMAG_DIR: $REMAG_DIR"
log "BUSCO_DIR: $BUSCO_DIR"
log "Minimum bin size: ${MIN_SIZE} bp"

# ============================================================
# TOOL CHECK
# ============================================================

command -v busco >/dev/null 2>&1 || {
    log "ERROR: BUSCO not found"
    log "Make sure the BUSCO environment/module is loaded"
    exit 1
}

# ============================================================
# FUNCTION: RUN BUSCO FOR ONE BIN
# ============================================================

run_busco() {

    local BIN="$1"

    local SAMPLE
    local BIN_NAME
    local OUT_PREFIX
    local OUTDIR
    local SIZE

    SAMPLE=$(basename "$(dirname "$(dirname "$BIN")")")
    BIN_NAME=$(basename "$BIN" .fa)

    OUT_PREFIX="${SAMPLE}_${BIN_NAME}"
    OUTDIR="$BUSCO_DIR/$OUT_PREFIX"

    log "----------------------------------------"
    log "Sample: $SAMPLE"
    log "Bin: $BIN_NAME"
    log "Input: $BIN"
    log "Output: $OUTDIR"

    # --------------------------------------------------------
    # INPUT CHECK
    # --------------------------------------------------------

    [[ -s "$BIN" ]] || {
        log "WARNING: Bin missing or empty -> skipping"
        return 0
    }

    # --------------------------------------------------------
    # CALCULATE TRUE FASTA SEQUENCE LENGTH
    # --------------------------------------------------------
    # This counts nucleotide sequence length rather than
    # filesystem byte size.

    SIZE=$(awk '
        !/^>/ {
            gsub(/[[:space:]]/, "")
            n += length($0)
        }
        END {
            print n+0
        }
    ' "$BIN")

    log "Bin sequence length: ${SIZE} bp"

    # --------------------------------------------------------
    # FILTER SMALL BINS
    # --------------------------------------------------------

    if (( SIZE < MIN_SIZE )); then
        log "Bin smaller than ${MIN_SIZE} bp -> skipping BUSCO"
        return 0
    fi

    # --------------------------------------------------------
    # SKIP COMPLETED
    # --------------------------------------------------------

    if compgen -G "$OUTDIR/short_summary.*.txt" > /dev/null \
       && [[ "${FORCE:-0}" != "1" ]]; then

        log "BUSCO already completed -> skipping"
        return 0

    fi

    # --------------------------------------------------------
    # CLEAN OLD / PARTIAL RUN
    # --------------------------------------------------------

    if [[ -d "$OUTDIR" ]]; then
        log "Removing previous BUSCO directory:"
        log "$OUTDIR"
        rm -rf "$OUTDIR"
    fi

    # --------------------------------------------------------
    # RUN BUSCO
    # --------------------------------------------------------

    log "Running BUSCO"

    busco \
        -i "$BIN" \
        -l eukaryota_odb10 \
        -o "$OUT_PREFIX" \
        -m genome \
        -c "$THREADS" \
        --out_path "$BUSCO_DIR" \
        &>> "$LOG_DIR/${OUT_PREFIX}.log"

    # --------------------------------------------------------
    # VALIDATE
    # --------------------------------------------------------

    if ! compgen -G "$OUTDIR/short_summary.*.txt" > /dev/null; then

        log "ERROR: BUSCO failed for:"
        log "$OUT_PREFIX"
        log "Check:"
        log "$LOG_DIR/${OUT_PREFIX}.log"

        return 1
    fi

    log "BUSCO completed: $OUT_PREFIX"
}


# ============================================================
# SAMPLE / BIN SELECTION
# ============================================================

if [[ -n "${BIN_FILE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-BIN MODE
    # --------------------------------------------------------

    log "Running in single-bin mode"

    [[ -s "$BIN_FILE" ]] || {
        log "ERROR: BIN_FILE not found or empty:"
        log "$BIN_FILE"
        exit 1
    }

    run_busco "$BIN_FILE"


elif [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------
    # Run all reMAG bins belonging to this sample.

    log "Running in single-sample mode: $SAMPLE"

    SAMPLE_BIN_DIR="$REMAG_DIR/$SAMPLE/bins"

    [[ -d "$SAMPLE_BIN_DIR" ]] || {
        log "WARNING: No bins directory for $SAMPLE"
        log "$SAMPLE_BIN_DIR"
        exit 0
    }

    mapfile -t BIN_FILES < <(
        find "$SAMPLE_BIN_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.fa" \
            -size +0c \
            | sort
    )

    N=${#BIN_FILES[@]}

    log "Bins detected for $SAMPLE: $N"

    if [[ "$N" -eq 0 ]]; then
        log "No bins found for $SAMPLE -> nothing to run"
        exit 0
    fi

    for BIN in "${BIN_FILES[@]}"; do
        run_busco "$BIN"
    done


else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t BIN_FILES < <(
        find "$REMAG_DIR" \
            -mindepth 3 \
            -maxdepth 3 \
            -type f \
            -path "*/bins/*.fa" \
            -size +0c \
            | sort
    )

    N=${#BIN_FILES[@]}

    log "Total bins detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No reMAG bins found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job, export SAMPLE=<sample>, or set BIN_FILE=<file>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT BIN --------------------

    BIN="${BIN_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

    log "Array task: ${SLURM_ARRAY_TASK_ID} / $N"

    run_busco "$BIN"

fi

log "BUSCO step completed"

;;

# ============================================================
### BUSCO CLASSIFICATION (EUK MAG QUALITY)
# ============================================================
busco_classify)

set -euo pipefail
shopt -s nullglob

log "Starting BUSCO classification step"

BUSCO_DIR="$OUTPUT_DIR/busco"
REMAG_DIR="$OUTPUT_DIR/remag"

CLASS_DIR="$OUTPUT_DIR/busco_classified"
SUMMARY="$CLASS_DIR/busco_summary.tsv"

mkdir -p "$CLASS_DIR"/{HQ,MQ,LQ}

# -------------------- INIT SUMMARY --------------------
echo -e "Sample\tBin\tCompleteness\tSingle\tDuplicated\tFragmented\tMissing\tCategory" > "$SUMMARY"

# -------------------- COLLECT BUSCO FILES --------------------
mapfile -t BUSCO_FILES < <(
    for f in "$BUSCO_DIR"/*/short_summary.*.txt; do
        [[ -s "$f" ]] && echo "$f"
    done
)

N=${#BUSCO_FILES[@]}
log "Detected BUSCO summaries: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No BUSCO results found"
    exit 1
}

# -------------------- PARSE --------------------
for FILE in "${BUSCO_FILES[@]}"; do

    DIR=$(dirname "$FILE")
    PREFIX=$(basename "$DIR")

    SAMPLE="${PREFIX%%_*}"
    BIN_NAME="${PREFIX#*_}"

    BIN="$REMAG_DIR/$SAMPLE/bins/${BIN_NAME}.fa"

    # -------------------- EXTRACT METRICS --------------------
    LINE=$(grep "C:" "$FILE" || true)

    [[ -n "$LINE" ]] || {
        log "WARNING: Could not parse BUSCO line for $PREFIX"
        continue
    }

    C=$(echo "$LINE" | sed -E 's/.*C:([0-9.]+)%.*/\1/')
    S=$(echo "$LINE" | sed -E 's/.*S:([0-9.]+)%.*/\1/')
    D=$(echo "$LINE" | sed -E 's/.*D:([0-9.]+)%.*/\1/')
    F=$(echo "$LINE" | sed -E 's/.*F:([0-9.]+)%.*/\1/')
    M=$(echo "$LINE" | sed -E 's/.*M:([0-9.]+)%.*/\1/')

    # -------------------- CLASSIFICATION (CORRECT) --------------------
    CATEGORY=$(awk -v c="$C" -v d="$D" 'BEGIN {
        if (c >= 90 && d <= 5) print "HQ";
        else if (c >= 50 && d <= 10) print "MQ";
        else print "LQ";
    }')

    # -------------------- WRITE SUMMARY --------------------
    echo -e "${SAMPLE}\t${BIN_NAME}\t${C}\t${S}\t${D}\t${F}\t${M}\t${CATEGORY}" >> "$SUMMARY"

    # -------------------- COPY BIN --------------------
    if [[ -s "$BIN" ]]; then
        cp "$BIN" "$CLASS_DIR/$CATEGORY/${SAMPLE}_${BIN_NAME}.fa"
    else
        log "WARNING: Bin not found for $PREFIX"
    fi

done

# -------------------- SUMMARY COUNTS --------------------
HQ=$(find "$CLASS_DIR/HQ" -name "*.fa" 2>/dev/null | wc -l)
MQ=$(find "$CLASS_DIR/MQ" -name "*.fa" 2>/dev/null | wc -l)
LQ=$(find "$CLASS_DIR/LQ" -name "*.fa" 2>/dev/null | wc -l)

log "Classification summary:"
log "HQ bins: $HQ"
log "MQ bins: $MQ"
log "LQ bins: $LQ"

log "Summary table: $SUMMARY"

;;

# ============================================================
# BUSCO SUMMARY PER SAMPLE (AGGREGATED)
# ============================================================
busco_summary_sample)

set -euo pipefail

log() { echo "[$(date '+%F %T')] $*"; }

BUSCO="$OUTPUT_DIR/busco_classified/busco_summary.tsv"
OUT="$OUTPUT_DIR/busco_classified/busco_summary_per_sample.tsv"

log "Creating BUSCO per-sample summary"

awk '
BEGIN {
    FS="\t"; OFS="\t"
}

NR==1 { next }

{
    sample = $1

    comp[sample] += $3
    frag[sample] += $6
    miss[sample] += $7

    n[sample]++
}

END {
    print "sample","mean_completeness","mean_fragmented","mean_missing","n_bins"

    for (s in comp) {
        print s,
              comp[s]/n[s],
              frag[s]/n[s],
              miss[s]/n[s],
              n[s]
    }
}
' "$BUSCO" > "$OUT"

log "Created: $OUT"

;;

# ============================================================
# MERGE QUAST + BUSCO (PER SAMPLE)
# ============================================================
quast_busco_merged)

set -euo pipefail

log() { echo "[$(date '+%F %T')] $*"; }

QUAST="$OUTPUT_DIR/quast/quast_summary.tsv"
BUSCO="$OUTPUT_DIR/busco_classified/busco_summary_per_sample.tsv"
OUT="$OUTPUT_DIR/quast/quast_busco_merged.tsv"

log "Merging QUAST + BUSCO (final robust)"

# -------------------- PREP FILES --------------------
# remove headers + sort by sample
tail -n +2 "$QUAST" | sort -k1,1 > tmp.quast
tail -n +2 "$BUSCO" | sort -k1,1 > tmp.busco

# -------------------- JOIN + CLEAN --------------------
join -t $'\t' -a1 -e "NA" -o auto tmp.quast tmp.busco \
| awk '
BEGIN {
    FS="\t"; OFS="\t"

    print "sample","sample_code","incubation",
          "n_contigs","total_length","N50","GC",
          "mean_completeness","mean_fragmented","mean_missing","n_bins"
}

{
    # -------------------- QUAST --------------------
    sample     = $1
    code       = $2
    incub      = $3
    contigs    = $4
    totlen     = $5
    N50        = $6
    GC         = $7

    # -------------------- BUSCO (SAFE NUMERIC) --------------------
    comp = ($8  == "NA" || $8  == "" ? 0 : $8  + 0)
    frag = ($9  == "NA" || $9  == "" ? 0 : $9  + 0)
    miss = ($10 == "NA" || $10 == "" ? 0 : $10 + 0)
    bins = ($11 == "NA" || $11 == "" ? 0 : $11 + 0)

    # -------------------- PRINT CLEAN --------------------
    printf "%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.6f\t%.6f\t%.6f\t%d\n",
           sample, code, incub,
           contigs, totlen, N50, GC,
           comp, frag, miss, bins
}
' > "$OUT"

# -------------------- CLEAN --------------------
rm -f tmp.quast tmp.busco

log "Created: $OUT"

;;


# ============================================================
### CAT (CONTIG-LEVEL — TIARA EUKARYOTIC CONTIGS)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
cat_contigs)

set -euo pipefail
shopt -s nullglob

log "Starting CAT (contig-level, eukaryotic contigs)"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-16}"

TIARA_DIR="$OUTPUT_DIR/tiara"
CAT_DIR="$OUTPUT_DIR/cat_contigs"
LOG_DIR="$OUTPUT_DIR/logs/cat_contigs"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$CAT_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "TIARA_DIR: $TIARA_DIR"
log "CAT_DIR: $CAT_DIR"

# ============================================================
# INPUT DIRECTORY CHECK
# ============================================================

[[ -d "$TIARA_DIR" ]] || {
    log "ERROR: Tiara directory not found:"
    log "$TIARA_DIR"
    exit 1
}

# ============================================================
# DEPENDENCY CHECK
# ============================================================

command -v "$CAT_CMD" >/dev/null 2>&1 || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    log "ERROR: DIAMOND not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_cat() {

    local SAMPLE="$1"

    local SAMPLE_DIR="$TIARA_DIR/$SAMPLE"
    local CONTIGS="$SAMPLE_DIR/euk.fa"

    local OUTDIR="$CAT_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}"
    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local EXIT_CODE

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "Eukaryotic contigs: $CONTIGS"
    log "Output prefix: $PREFIX"

    # ========================================================
    # INPUT CHECK
    # ========================================================

    if [[ ! -d "$SAMPLE_DIR" ]]; then
        log "ERROR: Tiara sample directory not found:"
        log "$SAMPLE_DIR"
        return 1
    fi

    # Empty euk.fa can be a legitimate Tiara result
    if [[ ! -e "$CONTIGS" ]]; then
        log "ERROR: euk.fa does not exist:"
        log "$CONTIGS"
        return 1
    fi

    if [[ ! -s "$CONTIGS" ]]; then
        log "No eukaryotic contigs for $SAMPLE -> skipping CAT"
        return 0
    fi

    # ========================================================
    # SKIP IF DONE
    # ========================================================

    if [[ -s "${PREFIX}.contig2classification.txt" &&
          "${FORCE:-0}" != "1" ]]; then

        log "CAT already completed -> skipping"
        log "${PREFIX}.contig2classification.txt"

        return 0
    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then
        log "Removing previous CAT directory:"
        log "$OUTDIR"
        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    rm -f "$LOGFILE"

    # ========================================================
    # RUN CAT
    # ========================================================

    log "Running CAT_pack contigs"

    set +e

    "$CAT_CMD" contigs \
        -c "$CONTIGS" \
        -d "$CAT_DB" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --nproc "$THREADS" \
        --block_size 2 \
        2>&1 | tee "$LOGFILE"

    EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "CAT exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if [[ "$EXIT_CODE" -ne 0 ]]; then

        log "ERROR: CAT failed for $SAMPLE"
        log "Last lines of CAT log:"

        tail -n 50 "$LOGFILE" || true

        return 1
    fi

    # ========================================================
    # VALIDATE
    # ========================================================

    if [[ ! -s "${PREFIX}.contig2classification.txt" ]]; then

        log "ERROR: CAT completed but classification file is missing:"
        log "${PREFIX}.contig2classification.txt"

        return 1
    fi

    log "Finished CAT for sample: $SAMPLE"
    log "Classification:"
    log "${PREFIX}.contig2classification.txt"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility:
# if TARGET_SAMPLE is provided but SAMPLE is not, use it.
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"

    run_cat "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t SAMPLE_DIRS < <(
        find "$TIARA_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    log "Total Tiara sample directories detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No Tiara sample directories found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------

    SAMPLE_DIR="${SAMPLE_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$SAMPLE_DIR")

    log "Array task: ${SLURM_ARRAY_TASK_ID} / $N"

    run_cat "$SAMPLE"

fi

log "CAT eukaryotic contig step completed"

;;

# ============================================================
### CONTIG LENGTH EXTRACTION (TIARA EUK CONTIGS)
# ============================================================
contig_lengths)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting contig length extraction"

THREADS="${THREADS:-4}"
TARGET_SAMPLE="${TARGET_SAMPLE:-}"
FORCE="${FORCE:-0}"

TIARA_DIR="$OUTPUT_DIR/tiara"
LEN_DIR="$OUTPUT_DIR/contig_lengths"
LOG_DIR="$OUTPUT_DIR/logs/contig_lengths"

mkdir -p "$LEN_DIR" "$LOG_DIR"

# ============================================================
# LOAD CONDA ENVIRONMENT (SEQKIT)
# ============================================================

log "Loading conda environment (seqkit)"

# Initialize conda (important on HPC)
source /staging/leuven/stg_00100/Rodgee/miniforge3/etc/profile.d/conda.sh

# Activate environment
conda activate seqkit_env

# Check
command -v seqkit >/dev/null || { log "ERROR: seqkit not found in conda env"; exit 1; }

log "Using seqkit: $(which seqkit)"

# -------------------- DEPENDENCY CHECK --------------------
command -v seqkit >/dev/null || { log "ERROR: seqkit not found"; exit 1; }

# ============================================================
# FUNCTION TO PROCESS ONE SAMPLE
# ============================================================
run_lengths() {

    local SAMPLE="$1"
    local COND="$TIARA_DIR/$SAMPLE"

    CONTIGS="$COND/euk.fa"
    OUT="$LEN_DIR/${SAMPLE}_lengths.tsv"

    log "----------------------------------------"
    log "Processing: $SAMPLE"

    [[ -s "$CONTIGS" ]] || {
        log "No euk.fa → skip"
        return
    }

    # -------------------- SKIP IF DONE --------------------
    if [[ -s "$OUT" && "$FORCE" != "1" ]]; then
        log "Already done → skipping"
        return
    fi

    # -------------------- EXTRACT LENGTHS --------------------
    log "Extracting contig lengths"

    (
        echo -e "contig\tlength"
        seqkit fx2tab -n -l "$CONTIGS"
    ) > "$OUT.tmp"

    # -------------------- SAFE WRITE --------------------
    if [[ -s "$OUT.tmp" ]]; then
        mv "$OUT.tmp" "$OUT"
        log "Created: $OUT"
    else
        log "ERROR: Failed to generate length file"
        rm -f "$OUT.tmp"
        exit 1
    fi
}

# ============================================================
# MODE SELECTION
# ============================================================

# -------------------- MANUAL MODE --------------------
if [[ -n "$TARGET_SAMPLE" ]]; then
    log "Manual mode → $TARGET_SAMPLE"
    run_lengths "$TARGET_SAMPLE"

# -------------------- ARRAY MODE --------------------
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    mapfile -t SAMPLES < <(
        for d in "$TIARA_DIR"/*; do
            [[ -d "$d" ]] || continue
            [[ -s "$d/euk.fa" ]] && basename "$d"
        done | sort
    )

    N=${#SAMPLES[@]}

    if (( SLURM_ARRAY_TASK_ID > N )); then
        log "No sample → exiting"
        exit 0
    fi

    SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
    run_lengths "$SAMPLE"

# -------------------- AUTO MODE --------------------
else
    log "Auto mode → processing all samples"

    for d in "$TIARA_DIR"/*; do
        [[ -d "$d" ]] || continue

        SAMPLE=$(basename "$d")
        run_lengths "$SAMPLE"
    done
fi

log "Contig length extraction completed"
;;

# ============================================================
### CAT (CONTIG-LEVEL — TIARA MITOCHONDRIAL CONTIGS)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
cat_mito_contigs)

set -euo pipefail
shopt -s nullglob

log "Starting CAT (mitochondrial contigs) step"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-16}"

TIARA_DIR="$OUTPUT_DIR/tiara"
CAT_DIR="$OUTPUT_DIR/cat_mito"
LOG_DIR="$OUTPUT_DIR/logs/cat_mito"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$CAT_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "TIARA_DIR: $TIARA_DIR"
log "CAT_DIR: $CAT_DIR"

# ============================================================
# INPUT DIRECTORY CHECK
# ============================================================

[[ -d "$TIARA_DIR" ]] || {
    log "ERROR: Tiara directory not found:"
    log "$TIARA_DIR"
    exit 1
}

# ============================================================
# DEPENDENCY CHECK
# ============================================================

command -v "$CAT_CMD" >/dev/null 2>&1 || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    log "ERROR: DIAMOND not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_cat_mito() {

    local SAMPLE="$1"

    local SAMPLE_DIR="$TIARA_DIR/$SAMPLE"
    local CONTIGS="$SAMPLE_DIR/mito.fa"

    local OUTDIR="$CAT_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}_mito"
    local LOGFILE="$LOG_DIR/${SAMPLE}_mito.log"

    local EXIT_CODE

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "Mitochondrial contigs: $CONTIGS"
    log "Output prefix: $PREFIX"

    # ========================================================
    # INPUT CHECK
    # ========================================================

    if [[ ! -d "$SAMPLE_DIR" ]]; then
        log "ERROR: Tiara sample directory not found:"
        log "$SAMPLE_DIR"
        return 1
    fi

    # Missing mito.fa indicates an incomplete Tiara result
    if [[ ! -e "$CONTIGS" ]]; then
        log "ERROR: mito.fa does not exist:"
        log "$CONTIGS"
        return 1
    fi

    # Empty mito.fa is a valid Tiara result
    if [[ ! -s "$CONTIGS" ]]; then
        log "No mitochondrial contigs for $SAMPLE -> skipping CAT"
        return 0
    fi

    # ========================================================
    # SKIP IF DONE
    # ========================================================

    if [[ -s "${PREFIX}.contig2classification.txt" &&
          "${FORCE:-0}" != "1" ]]; then

        log "CAT mitochondrial classification already completed -> skipping"
        log "${PREFIX}.contig2classification.txt"

        return 0
    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then
        log "Removing previous CAT mitochondrial directory:"
        log "$OUTDIR"
        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    rm -f "$LOGFILE"

    # ========================================================
    # RUN CAT
    # ========================================================

    log "Running CAT_pack on mitochondrial contigs"

    set +e

    "$CAT_CMD" contigs \
        -c "$CONTIGS" \
        -d "$CAT_DB" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --nproc "$THREADS" \
        --block_size 2 \
        2>&1 | tee "$LOGFILE"

    EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "CAT exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if [[ "$EXIT_CODE" -ne 0 ]]; then

        log "ERROR: CAT failed for mitochondrial contigs of $SAMPLE"
        log "Last lines of CAT log:"

        tail -n 50 "$LOGFILE" || true

        return 1
    fi

    # ========================================================
    # VALIDATE
    # ========================================================

    if [[ ! -s "${PREFIX}.contig2classification.txt" ]]; then

        log "ERROR: CAT completed but classification file is missing:"
        log "${PREFIX}.contig2classification.txt"

        return 1
    fi

    log "Finished CAT mitochondrial classification for $SAMPLE"
    log "Classification:"
    log "${PREFIX}.contig2classification.txt"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"

    run_cat_mito "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t SAMPLE_DIRS < <(
        find "$TIARA_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    log "Total Tiara sample directories detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No Tiara sample directories found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------

    SAMPLE_DIR="${SAMPLE_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$SAMPLE_DIR")

    log "Array task: ${SLURM_ARRAY_TASK_ID} / $N"

    run_cat_mito "$SAMPLE"

fi

log "CAT mitochondrial contig step completed"

;;

# ============================================================
### BAT (BIN-LEVEL — REMAG EUKARYOTIC BINS)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
bat_bins)

set -euo pipefail
shopt -s nullglob

log "Starting BAT (bin-level, reMAG eukaryotic bins)"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-16}"

REMAG_DIR="$OUTPUT_DIR/remag"
BAT_DIR="$OUTPUT_DIR/bat_bins"
LOG_DIR="$OUTPUT_DIR/logs/bat_bins"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$BAT_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "REMAG_DIR: $REMAG_DIR"
log "BAT_DIR: $BAT_DIR"

# ============================================================
# INPUT DIRECTORY CHECK
# ============================================================

[[ -d "$REMAG_DIR" ]] || {
    log "ERROR: reMAG directory not found:"
    log "$REMAG_DIR"
    exit 1
}

# ============================================================
# DEPENDENCY CHECK
# ============================================================

command -v "$CAT_CMD" >/dev/null 2>&1 || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    log "ERROR: DIAMOND not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_bat() {

    local SAMPLE="$1"

    local SAMPLE_DIR="$REMAG_DIR/$SAMPLE"
    local BINS_DIR="$SAMPLE_DIR/bins"

    local OUTDIR="$BAT_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}"
    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local EXIT_CODE
    local NBINS

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "reMAG directory: $SAMPLE_DIR"
    log "Bins directory: $BINS_DIR"
    log "Output prefix: $PREFIX"

    # ========================================================
    # INPUT CHECK
    # ========================================================

    if [[ ! -d "$SAMPLE_DIR" ]]; then
        log "ERROR: reMAG sample directory not found:"
        log "$SAMPLE_DIR"
        return 1
    fi

    # No bins is a valid reMAG result
    if [[ ! -d "$BINS_DIR" ]]; then
        log "No bins directory for $SAMPLE -> skipping BAT"
        return 0
    fi

    mapfile -t BIN_FILES < <(
        find "$BINS_DIR" \
            -maxdepth 1 \
            -type f \
            -name "*.fa" \
            -size +0c \
            | sort
    )

    NBINS=${#BIN_FILES[@]}

    if [[ "$NBINS" -eq 0 ]]; then
        log "No non-empty reMAG bins for $SAMPLE -> skipping BAT"
        return 0
    fi

    log "Number of bins: $NBINS"

    # ========================================================
    # SKIP IF DONE
    # ========================================================

    if [[ -s "${PREFIX}.bin2classification.txt" &&
          "${FORCE:-0}" != "1" ]]; then

        log "BAT already completed -> skipping"
        log "${PREFIX}.bin2classification.txt"

        return 0
    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then
        log "Removing previous BAT directory:"
        log "$OUTDIR"
        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    rm -f "$LOGFILE"

    # ========================================================
    # RUN BAT
    # ========================================================

    log "Running CAT_pack bins"

    set +e

    "$CAT_CMD" bins \
        -b "$BINS_DIR" \
        -d "$CAT_DB" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --bin_suffix .fa \
        --nproc "$THREADS" \
        --block_size 2 \
        2>&1 | tee "$LOGFILE"

    EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "BAT exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if [[ "$EXIT_CODE" -ne 0 ]]; then

        log "ERROR: BAT failed for $SAMPLE"
        log "Last lines of BAT log:"

        tail -n 50 "$LOGFILE" || true

        return 1
    fi

    # ========================================================
    # VALIDATE
    # ========================================================

    if [[ ! -s "${PREFIX}.bin2classification.txt" ]]; then

        log "ERROR: BAT completed but classification file is missing:"
        log "${PREFIX}.bin2classification.txt"

        return 1
    fi

    log "Finished BAT for sample: $SAMPLE"
    log "Classification:"
    log "${PREFIX}.bin2classification.txt"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"

    run_bat "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t SAMPLE_DIRS < <(
        find "$REMAG_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    log "Total reMAG sample directories detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No reMAG sample directories found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------

    SAMPLE_DIR="${SAMPLE_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$SAMPLE_DIR")

    log "Array task: ${SLURM_ARRAY_TASK_ID} / $N"

    run_bat "$SAMPLE"

fi

log "BAT bin-level step completed"

;;

# ============================================================
### RAT (READ-LEVEL — TIARA EUK CONTIGS + CAT OUTPUTS)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
rat_reads)

set -euo pipefail
shopt -s nullglob

log "Starting RAT (read-level classification)"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-8}"

export OMP_NUM_THREADS="$THREADS"

TRIM_DIR="$OUTPUT_DIR/trimmed"
TIARA_DIR="$OUTPUT_DIR/tiara"
TIARA_CLEAN_DIR="$OUTPUT_DIR/tiara_clean"

CAT_CONTIG_DIR="$OUTPUT_DIR/cat_contigs"

RAT_DIR="$OUTPUT_DIR/rat_reads"
LOG_DIR="$OUTPUT_DIR/logs/rat_reads"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p \
    "$RAT_DIR" \
    "$LOG_DIR" \
    "$TIARA_CLEAN_DIR"

log "THREADS: $THREADS"
log "TRIM_DIR: $TRIM_DIR"
log "TIARA_DIR: $TIARA_DIR"
log "CAT_CONTIG_DIR: $CAT_CONTIG_DIR"
log "RAT_DIR: $RAT_DIR"

# ============================================================
# INPUT DIRECTORY CHECKS
# ============================================================

[[ -d "$TRIM_DIR" ]] || {
    log "ERROR: Trimmed-read directory not found:"
    log "$TRIM_DIR"
    exit 1
}

[[ -d "$TIARA_DIR" ]] || {
    log "ERROR: Tiara directory not found:"
    log "$TIARA_DIR"
    exit 1
}

[[ -d "$CAT_CONTIG_DIR" ]] || {
    log "ERROR: CAT contig directory not found:"
    log "$CAT_CONTIG_DIR"
    exit 1
}

[[ -d "$CAT_TAX" ]] || {
    log "ERROR: CAT taxonomy directory not found:"
    log "$CAT_TAX"
    exit 1
}

# ============================================================
# DEPENDENCY CHECK
# ============================================================

command -v "$CAT_CMD" >/dev/null 2>&1 || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    log "ERROR: DIAMOND not found"
    exit 1
}

command -v bwa >/dev/null 2>&1 || {
    log "ERROR: BWA not found"
    exit 1
}

command -v samtools >/dev/null 2>&1 || {
    log "ERROR: samtools not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_rat() {

    local SAMPLE="$1"

    local R1="$TRIM_DIR/${SAMPLE}_trimmed_R1.fastq.gz"
    local R2="$TRIM_DIR/${SAMPLE}_trimmed_R2.fastq.gz"

    local RAW_CONTIGS="$TIARA_DIR/$SAMPLE/euk.fa"

    local CLEAN_DIR="$TIARA_CLEAN_DIR/$SAMPLE"
    local CONTIGS="$CLEAN_DIR/euk_clean.fa"

    local CAT_SAMPLE_DIR="$CAT_CONTIG_DIR/$SAMPLE"

    local RAW_C2C
    local CAT_C2C="$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.clean.txt"

    local OUTDIR="$RAT_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}"

    local RAT_TMPDIR="$OUTDIR/tmp"
    local DUMMY_U2C="$OUTDIR/dummy_u2c.txt"

    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local COMPLETE_FILE="${PREFIX}.RAT.complete.abundance.txt"
    local CONTIG_FILE="${PREFIX}.RAT.contig.abundance.txt"
    local READ_FILE="${PREFIX}.RAT.read2classification.txt"

    local EXIT_CODE

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "R1: $R1"
    log "R2: $R2"
    log "Raw euk contigs: $RAW_CONTIGS"
    log "RAT output: $OUTDIR"

    # ========================================================
    # INPUT CHECK: READS
    # ========================================================

    [[ -s "$R1" ]] || {
        log "ERROR: Missing R1:"
        log "$R1"
        return 1
    }

    [[ -s "$R2" ]] || {
        log "ERROR: Missing R2:"
        log "$R2"
        return 1
    }

    # ========================================================
    # INPUT CHECK: TIARA EUKARYOTIC CONTIGS
    # ========================================================

    if [[ ! -e "$RAW_CONTIGS" ]]; then
        log "ERROR: Tiara euk.fa does not exist:"
        log "$RAW_CONTIGS"
        return 1
    fi

    # Empty euk.fa is a legitimate Tiara result
    if [[ ! -s "$RAW_CONTIGS" ]]; then
        log "No eukaryotic contigs for $SAMPLE -> skipping RAT"
        return 0
    fi

    # ========================================================
    # CAT CONTIG CLASSIFICATION
    # ========================================================
    #
    # Support both possible CAT filename conventions:
    #   sample.contig2classification.txt
    #   sample.contig2classification.txt
    #

    if [[ -s "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt" ]]; then

        RAW_C2C="$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"

    elif [[ -s "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt" ]]; then

        RAW_C2C="$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"

    else

        log "ERROR: CAT contig classification not found for $SAMPLE"
        log "Checked:"
        log "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"
        log "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"

        return 1
    fi

    log "CAT classification: $RAW_C2C"

    # ========================================================
    # SKIP IF COMPLETED
    # ========================================================

    if [[ -s "$COMPLETE_FILE" &&
          -s "$CONTIG_FILE" &&
          -s "$READ_FILE" &&
          "${FORCE:-0}" != "1" ]]; then

        log "RAT already completed -> skipping"
        return 0
    fi

    # ========================================================
    # CLEAN OLD / PARTIAL RAT OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then
        log "Removing previous RAT output:"
        log "$OUTDIR"
        rm -rf "$OUTDIR"
    fi

    mkdir -p \
        "$OUTDIR" \
        "$RAT_TMPDIR" \
        "$CLEAN_DIR"

    rm -f "$LOGFILE"

    # ========================================================
    # CLEAN FASTA HEADERS
    # ========================================================

    log "Cleaning FASTA headers"

    awk '
    /^>/ {
        sub(/^>/, "")
        split($0, a, " ")
        print ">" a[1]
        next
    }
    {
        print
    }
    ' "$RAW_CONTIGS" > "$CONTIGS"

    [[ -s "$CONTIGS" ]] || {
        log "ERROR: Cleaned contig FASTA is empty"
        return 1
    }

    # ========================================================
    # CLEAN CAT CONTIG CLASSIFICATION
    # ========================================================

    log "Cleaning CAT contig classification"

    awk '
    BEGIN {
        FS=OFS="\t"
    }

    /^#/ {
        print
        next
    }

    {
        split($1, a, " ")
        $1=a[1]
        print
    }
    ' "$RAW_C2C" > "$CAT_C2C"

    [[ -s "$CAT_C2C" ]] || {
        log "ERROR: Cleaned CAT classification is empty"
        return 1
    }

    # ========================================================
    # DUMMY U2C FILE
    # ========================================================

    touch "$DUMMY_U2C"

    # ========================================================
    # RUN RAT
    # ========================================================

    log "Running CAT_pack reads (RAT)"

    set +e

    "$CAT_CMD" reads \
        --mode cr \
        -c "$CONTIGS" \
        --c2c "$CAT_C2C" \
        --u2c "$DUMMY_U2C" \
        -1 "$R1" \
        -2 "$R2" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --tmpdir "$RAT_TMPDIR" \
        --read2classification \
        --force \
        --nproc "$THREADS" \
        2>&1 | tee "$LOGFILE"

    EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "RAT exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if [[ "$EXIT_CODE" -ne 0 ]]; then

        log "ERROR: RAT failed for $SAMPLE"
        log "Last lines of RAT log:"

        tail -n 50 "$LOGFILE" || true

        rm -f "$DUMMY_U2C"

        return 1
    fi

    # ========================================================
    # VALIDATE OUTPUT
    # ========================================================

    [[ -s "$READ_FILE" ]] || {
        log "ERROR: Missing RAT read2classification output"
        rm -f "$DUMMY_U2C"
        return 1
    }

    [[ -s "$COMPLETE_FILE" ]] || {
        log "ERROR: Missing RAT complete abundance output"
        rm -f "$DUMMY_U2C"
        return 1
    }

    [[ -s "$CONTIG_FILE" ]] || {
        log "ERROR: Missing RAT contig abundance output"
        rm -f "$DUMMY_U2C"
        return 1
    }

    # ========================================================
    # SUMMARY
    # ========================================================

    log "RAT outputs generated:"

    ls -lh "$OUTDIR"/*.RAT.* 2>/dev/null || true

    # ========================================================
    # CLEANUP
    # ========================================================

    rm -f "$DUMMY_U2C"

    if [[ "${DEBUG:-0}" != "1" ]]; then
        rm -rf "$RAT_TMPDIR"
    else
        log "DEBUG=1 -> keeping RAT temporary directory:"
        log "$RAT_TMPDIR"
    fi

    log "Finished RAT for $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility:
# TARGET_SAMPLE still works, but SAMPLE is now preferred.
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"

    run_rat "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    # Use Tiara directories as the stable sample set.
    # Samples with an empty euk.fa will simply skip inside
    # run_rat().
    mapfile -t SAMPLE_DIRS < <(
        find "$TIARA_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    log "Total Tiara sample directories detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No Tiara sample directories found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------

    SAMPLE_DIR="${SAMPLE_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$SAMPLE_DIR")

    log "Array task: ${SLURM_ARRAY_TASK_ID} / $N"

    run_rat "$SAMPLE"

fi

log "RAT read-level step completed"

;;

# ============================================================
### EUKARYOTIC READ LENGTH EXTRACTION
### (READS MAPPING TO EUKARYOTIC CONTIGS)
# ============================================================
read_length_euk)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

RAT_DIR="$OUTPUT_DIR/rat_reads_summary"
MAP_DIR="$OUTPUT_DIR/mapping"
LEN_DIR="$OUTPUT_DIR/read_length_euk"
LOG_DIR="$OUTPUT_DIR/logs/read_length_euk"

mkdir -p "$LEN_DIR" "$LOG_DIR"

log "Starting eukaryotic read length extraction"

# ============================================================
# SELECT SAMPLE
# ============================================================

if [[ -n "${TARGET_SAMPLE:-}" ]]; then

    SAMPLE="$TARGET_SAMPLE"
    RAT_FILE="$RAT_DIR/${SAMPLE}.taxonomy_structured.tsv"

    log "Using TARGET_SAMPLE: $SAMPLE"

else

    mapfile -t RAT_FILES < <(find "$RAT_DIR" -name "*.taxonomy_structured.tsv" | sort)

    N=${#RAT_FILES[@]}

    log "Total samples: $N"

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID out of range (1-$N)"
        exit 1
    fi

    RAT_FILE="${RAT_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$RAT_FILE" .taxonomy_structured.tsv)

fi

# ============================================================
# SAMPLE-SPECIFIC FILES
# ============================================================

BAM="$MAP_DIR/${SAMPLE}.bam"

EUK_CONTIGS="$LEN_DIR/${SAMPLE}.euk.contigs.txt"
BAM_REFS="$LEN_DIR/${SAMPLE}.bam_refs.txt"
FULL_REFS="$LEN_DIR/${SAMPLE}.full_refs.txt"

TMP_DIR="$LEN_DIR/${SAMPLE}_tmp_refs"

EUK_BAM="$LEN_DIR/${SAMPLE}.euk.bam"
EUK_FASTQ="$LEN_DIR/${SAMPLE}.euk.fastq"
LEN_ALL="$LEN_DIR/${SAMPLE}.all_lengths.txt"
OUT="$LEN_DIR/${SAMPLE}.read_length.tsv"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# ============================================================
# INPUT CHECKS
# ============================================================

[[ -s "$RAT_FILE" ]] || {
    log "ERROR: Missing RAT file"
    exit 1
}

[[ -s "$BAM" ]] || {
    log "ERROR: Missing BAM file"
    exit 1
}

samtools quickcheck "$BAM"

# ============================================================
# SKIP
# ============================================================

if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
    log "Output exists -> skipping"
    exit 0
fi

# ============================================================
# STEP 1: EXTRACT EUKARYOTIC CONTIGS
# ============================================================

log "Extracting eukaryotic contigs"

awk '$5=="Eukaryota" {print $1}' \
    "$RAT_FILE" > "$EUK_CONTIGS"

N_CONTIGS=$(wc -l < "$EUK_CONTIGS")

log "Eukaryotic contigs: $N_CONTIGS"

[[ "$N_CONTIGS" -gt 0 ]] || {
    log "WARNING: No eukaryotic contigs found"
    exit 0
}

# ============================================================
# STEP 2: MATCH CONTIGS TO BAM REFERENCES
# ============================================================

log "Matching contigs to BAM references"

samtools idxstats "$BAM" | cut -f1 > "$BAM_REFS"

awk '
NR==FNR {
    keep[$1]=1
    next
}
{
    split($0,a," ")
    if(a[1] in keep)
        print $0
}
' "$EUK_CONTIGS" "$BAM_REFS" > "$FULL_REFS"

N_REFS=$(wc -l < "$FULL_REFS")

log "Matched BAM references: $N_REFS"

[[ "$N_REFS" -gt 0 ]] || {
    log "ERROR: No matching BAM references found"
    exit 1
}

# ============================================================
# STEP 3: EXTRACT EUKARYOTIC READS
# ============================================================

log "Extracting reads from eukaryotic contigs"

mapfile -t REFS < "$FULL_REFS"

log "Number of references: ${#REFS[@]}"

samtools view -b \
    "$BAM" \
    "${REFS[@]}" \
    > "$EUK_BAM"

[[ -s "$EUK_BAM" ]] || {
    log "ERROR: Euk BAM not created"
    exit 1
}

N_EUK=$(samtools view -c "$EUK_BAM")

log "Alignments retained: $N_EUK"

[[ "$N_EUK" -gt 0 ]] || {
    log "ERROR: No alignments retained"
    exit 1
}

# ============================================================
# STEP 4: BAM -> FASTQ
# ============================================================

log "Converting BAM to FASTQ"

samtools fastq "$EUK_BAM" > "$EUK_FASTQ"

N_FASTQ=$(grep -c '^@' "$EUK_FASTQ")

log "FASTQ reads: $N_FASTQ"

[[ "$N_FASTQ" -gt 0 ]] || {
    log "ERROR: FASTQ contains no reads"
    exit 1
}

# ============================================================
# STEP 5: READ LENGTH EXTRACTION
# ============================================================

log "Extracting read lengths"

awk 'NR%4==2 {print length($0)}' \
    "$EUK_FASTQ" > "$LEN_ALL"

N_LEN=$(wc -l < "$LEN_ALL")

log "Length entries: $N_LEN"

[[ "$N_LEN" -gt 0 ]] || {
    log "ERROR: No read lengths extracted"
    exit 1
}

# ============================================================
# STEP 6: SUMMARY STATISTICS
# ============================================================

log "Computing summary statistics"

awk -v sample="$SAMPLE" '
{
    sum += $1
    n++
    values[n] = $1
}
END{
    if(n==0){
        print sample "\t0\t0\t0"
        exit
    }

    mean = sum / n

    asort(values)

    if(n % 2)
        median = values[(n+1)/2]
    else
        median = (values[n/2] + values[n/2+1]) / 2

    print sample "\t" n "\t" mean "\t" median
}
' "$LEN_ALL" > "$OUT"

# ============================================================
# CLEANUP
# ============================================================

rm -f \
    "$BAM_REFS" \
    "$FULL_REFS" \
    "$EUK_FASTQ"

# Keep these for debugging / possible reruns
# rm -f "$EUK_BAM"
# rm -f "$LEN_ALL"

# ============================================================
# VALIDATE
# ============================================================

[[ -s "$OUT" ]] || {
    log "ERROR: Output not created"
    exit 1
}

log "Preview:"
cat "$OUT"

log "Finished euk read length extraction for $SAMPLE"
;;

# ============================================================
### CAT (CONTIG-LEVEL — TIARA PROKARYOTIC CONTIGS)
# ============================================================
cat_contigs_prok)

set -euo pipefail
shopt -s nullglob

log "Starting CAT (contig-level, prok)"

# ============================================================
# CONFIG
# ============================================================

THREADS="${THREADS:-8}"

TIARA_DIR="$OUTPUT_DIR/tiara"
CAT_DIR="$OUTPUT_DIR/cat_contigs_prok"
LOG_DIR="$OUTPUT_DIR/logs/cat_contigs_prok"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

CAT_CMD="${CAT_CMD:-CAT_pack}"

# Backward compatibility:
# if TARGET_SAMPLE was supplied but SAMPLE was not,
# use TARGET_SAMPLE as SAMPLE.
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

mkdir -p "$CAT_DIR" "$LOG_DIR"

log "THREADS: $THREADS"
log "TIARA_DIR: $TIARA_DIR"
log "CAT_DIR: $CAT_DIR"

# ============================================================
# INPUT DIRECTORY CHECK
# ============================================================

[[ -d "$TIARA_DIR" ]] || {
    log "ERROR: Tiara directory not found:"
    log "$TIARA_DIR"
    exit 1
}

# ============================================================
# DATABASE CHECKS
# ============================================================

[[ -d "$CAT_DB" ]] || {
    log "ERROR: CAT database directory not found:"
    log "$CAT_DB"
    exit 1
}

[[ -d "$CAT_TAX" ]] || {
    log "ERROR: CAT taxonomy directory not found:"
    log "$CAT_TAX"
    exit 1
}

# ============================================================
# DEPENDENCY CHECKS
# ============================================================

command -v "$CAT_CMD" >/dev/null 2>&1 || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    log "ERROR: DIAMOND not found"
    exit 1
}

log "CAT executable: $(command -v "$CAT_CMD")"
log "DIAMOND executable: $(command -v diamond)"

# ============================================================
# FUNCTION TO PROCESS ONE SAMPLE
# ============================================================

run_cat() {

    local SAMPLE="$1"

    local SAMPLE_DIR="$TIARA_DIR/$SAMPLE"
    local CONTIGS="$SAMPLE_DIR/prok.fa"

    local OUTDIR="$CAT_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}"

    local CAT_LOG="$LOG_DIR/${SAMPLE}.log"

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "Input: $CONTIGS"
    log "Output: $OUTDIR"

    # ========================================================
    # SAMPLE DIRECTORY CHECK
    # ========================================================

    [[ -d "$SAMPLE_DIR" ]] || {
        log "ERROR: Tiara sample directory not found:"
        log "$SAMPLE_DIR"
        return 1
    }

    # ========================================================
    # INPUT CHECK
    # ========================================================

    if [[ ! -e "$CONTIGS" ]]; then
        log "ERROR: prok.fa does not exist:"
        log "$CONTIGS"
        return 1
    fi

    if [[ ! -s "$CONTIGS" ]]; then
        log "No prokaryotic contigs detected for $SAMPLE"
        log "CAT classification not required -> skipping"
        return 0
    fi

    # ========================================================
    # SKIP COMPLETED
    # ========================================================

    if [[ -s "${PREFIX}.contig2classification.txt" &&
          "${FORCE:-0}" != "1" ]]; then

        log "CAT already completed -> skipping"
        log "${PREFIX}.contig2classification.txt"

        return 0
    fi

    # ========================================================
    # CLEAN OLD / PARTIAL OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then
        log "Removing previous CAT output:"
        log "$OUTDIR"
        rm -rf "$OUTDIR"
    fi

    mkdir -p "$OUTDIR"

    # Remove old sample log before rerunning
    rm -f "$CAT_LOG"

    # ========================================================
    # RUN CAT
    # ========================================================

    log "Running CAT_pack contigs"
    log "Threads: $THREADS"

    set +e

    "$CAT_CMD" contigs \
        -c "$CONTIGS" \
        -d "$CAT_DB" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --nproc "$THREADS" \
        --block_size 2 \
        2>&1 | tee "$CAT_LOG"

    local EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "CAT exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if [[ "$EXIT_CODE" -ne 0 ]]; then

        log "ERROR: CAT failed for $SAMPLE"
        log "Last 50 lines of CAT log:"

        tail -n 50 "$CAT_LOG" || true

        return "$EXIT_CODE"
    fi

    # ========================================================
    # VALIDATE OUTPUT
    # ========================================================

    if [[ ! -s "${PREFIX}.contig2classification.txt" ]]; then

        log "ERROR: CAT finished but classification output is missing:"
        log "${PREFIX}.contig2classification.txt"

        return 1
    fi

    log "Finished CAT for sample: $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"
    log "Sample: $SAMPLE"

    run_cat "$SAMPLE"

elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    mapfile -t SAMPLE_DIRS < <(
        find "$TIARA_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    log "Total Tiara sample directories detected: $N"

    if (( N == 0 )); then
        log "ERROR: No Tiara sample directories found"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Array task ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    SAMPLE_DIR="${SAMPLE_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$SAMPLE_DIR")

    log "Array task ${SLURM_ARRAY_TASK_ID} -> $SAMPLE"

    run_cat "$SAMPLE"

else

    # --------------------------------------------------------
    # LOCAL / SEQUENTIAL MODE
    # --------------------------------------------------------

    log "Running in local sequential mode"

    mapfile -t SAMPLE_DIRS < <(
        find "$TIARA_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    if (( N == 0 )); then
        log "ERROR: No Tiara sample directories found"
        exit 1
    fi

    log "Total Tiara samples detected: $N"

    for SAMPLE_DIR in "${SAMPLE_DIRS[@]}"; do

        SAMPLE=$(basename "$SAMPLE_DIR")

        run_cat "$SAMPLE"

    done

fi

log "CAT prok contig step completed"

;;

# ============================================================
### RAT (READ-LEVEL — TIARA PROK CONTIGS + CAT PROK OUTPUTS)
### ARRAY OR SINGLE-SAMPLE
# ============================================================
rat_reads_prok)

set -euo pipefail
shopt -s nullglob

log "Starting RAT (read-level classification, prokaryotic contigs)"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-8}"

export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

TRIM_DIR="$OUTPUT_DIR/trimmed"

TIARA_DIR="$OUTPUT_DIR/tiara"
TIARA_CLEAN_DIR="$OUTPUT_DIR/tiara_clean_prok"

CAT_CONTIG_DIR="$OUTPUT_DIR/cat_contigs_prok"

RAT_DIR="$OUTPUT_DIR/rat_reads_prok"
LOG_DIR="$OUTPUT_DIR/logs/rat_reads_prok"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p \
    "$RAT_DIR" \
    "$LOG_DIR" \
    "$TIARA_CLEAN_DIR"

log "THREADS: $THREADS"
log "TRIM_DIR: $TRIM_DIR"
log "TIARA_DIR: $TIARA_DIR"
log "CAT_CONTIG_DIR: $CAT_CONTIG_DIR"
log "RAT_DIR: $RAT_DIR"

# ============================================================
# INPUT DIRECTORY CHECKS
# ============================================================

[[ -d "$TRIM_DIR" ]] || {
    log "ERROR: Trimmed-read directory not found:"
    log "$TRIM_DIR"
    exit 1
}

[[ -d "$TIARA_DIR" ]] || {
    log "ERROR: Tiara directory not found:"
    log "$TIARA_DIR"
    exit 1
}

[[ -d "$CAT_CONTIG_DIR" ]] || {
    log "ERROR: CAT prok contig directory not found:"
    log "$CAT_CONTIG_DIR"
    exit 1
}

[[ -d "$CAT_TAX" ]] || {
    log "ERROR: CAT taxonomy directory not found:"
    log "$CAT_TAX"
    exit 1
}

# ============================================================
# DEPENDENCY CHECK
# ============================================================

command -v "$CAT_CMD" >/dev/null 2>&1 || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null 2>&1 || {
    log "ERROR: DIAMOND not found"
    exit 1
}

command -v bwa >/dev/null 2>&1 || {
    log "ERROR: BWA not found"
    exit 1
}

command -v samtools >/dev/null 2>&1 || {
    log "ERROR: samtools not found"
    exit 1
}

# ============================================================
# FUNCTION: PROCESS ONE SAMPLE
# ============================================================

run_rat_prok() {

    local SAMPLE="$1"

    local R1="$TRIM_DIR/${SAMPLE}_trimmed_R1.fastq.gz"
    local R2="$TRIM_DIR/${SAMPLE}_trimmed_R2.fastq.gz"

    local RAW_CONTIGS="$TIARA_DIR/$SAMPLE/prok.fa"

    local CLEAN_DIR="$TIARA_CLEAN_DIR/$SAMPLE"
    local CONTIGS="$CLEAN_DIR/prok_clean.fa"

    local CAT_SAMPLE_DIR="$CAT_CONTIG_DIR/$SAMPLE"

    local RAW_C2C
    local CAT_C2C="$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.clean.txt"

    local OUTDIR="$RAT_DIR/$SAMPLE"
    local PREFIX="$OUTDIR/${SAMPLE}"

    local RAT_TMPDIR="$OUTDIR/tmp"
    local DUMMY_U2C="$OUTDIR/dummy_u2c.txt"

    local LOGFILE="$LOG_DIR/${SAMPLE}.log"

    local COMPLETE_FILE="${PREFIX}.RAT.complete.abundance.txt"
    local CONTIG_FILE="${PREFIX}.RAT.contig.abundance.txt"
    local READ_FILE="${PREFIX}.RAT.read2classification.txt"

    local EXIT_CODE

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"
    log "R1: $R1"
    log "R2: $R2"
    log "Raw prok contigs: $RAW_CONTIGS"
    log "RAT output: $OUTDIR"

    # ========================================================
    # INPUT CHECK: READS
    # ========================================================

    [[ -s "$R1" ]] || {
        log "ERROR: Missing R1:"
        log "$R1"
        return 1
    }

    [[ -s "$R2" ]] || {
        log "ERROR: Missing R2:"
        log "$R2"
        return 1
    }

    # ========================================================
    # INPUT CHECK: TIARA PROKARYOTIC CONTIGS
    # ========================================================

    if [[ ! -e "$RAW_CONTIGS" ]]; then
        log "ERROR: Tiara prok.fa does not exist:"
        log "$RAW_CONTIGS"
        return 1
    fi

    # Empty prok.fa is a valid Tiara result
    if [[ ! -s "$RAW_CONTIGS" ]]; then
        log "No prokaryotic contigs for $SAMPLE -> skipping RAT"
        return 0
    fi

    # ========================================================
    # CAT PROK CONTIG CLASSIFICATION
    # ========================================================
    #
    # Support either CAT naming convention:
    #
    #   SAMPLE.contig2classification.txt
    #   SAMPLE.contig2classification.txt
    #

    if [[ -s "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt" ]]; then

        RAW_C2C="$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"

    elif [[ -s "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt" ]]; then

        RAW_C2C="$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"

    else

        log "ERROR: CAT prok contig classification not found for $SAMPLE"
        log "Checked:"
        log "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"
        log "$CAT_SAMPLE_DIR/${SAMPLE}.contig2classification.txt"

        return 1
    fi

    log "CAT classification: $RAW_C2C"

    # ========================================================
    # SKIP IF COMPLETED
    # ========================================================

    if [[ -s "$COMPLETE_FILE" &&
          -s "$CONTIG_FILE" &&
          -s "$READ_FILE" &&
          "${FORCE:-0}" != "1" ]]; then

        log "RAT prok already completed -> skipping"
        return 0
    fi

    # ========================================================
    # CLEAN OLD / PARTIAL RAT OUTPUT
    # ========================================================

    if [[ -d "$OUTDIR" ]]; then
        log "Removing previous RAT prok output:"
        log "$OUTDIR"
        rm -rf "$OUTDIR"
    fi

    mkdir -p \
        "$OUTDIR" \
        "$RAT_TMPDIR" \
        "$CLEAN_DIR"

    rm -f "$LOGFILE"

    # ========================================================
    # CLEAN FASTA HEADERS
    # ========================================================

    log "Cleaning prokaryotic FASTA headers"

    awk '
    /^>/ {
        sub(/^>/, "")
        split($0, a, " ")
        print ">" a[1]
        next
    }
    {
        print
    }
    ' "$RAW_CONTIGS" > "$CONTIGS"

    [[ -s "$CONTIGS" ]] || {
        log "ERROR: Cleaned prokaryotic FASTA is empty"
        return 1
    }

    # ========================================================
    # CLEAN CAT CONTIG CLASSIFICATION
    # ========================================================

    log "Cleaning CAT prok contig classification"

    awk '
    BEGIN {
        FS=OFS="\t"
    }

    /^#/ {
        print
        next
    }

    {
        split($1, a, " ")
        $1=a[1]
        print
    }
    ' "$RAW_C2C" > "$CAT_C2C"

    [[ -s "$CAT_C2C" ]] || {
        log "ERROR: Cleaned CAT classification is empty"
        return 1
    }

    # ========================================================
    # DUMMY U2C FILE
    # ========================================================

    touch "$DUMMY_U2C"

    # ========================================================
    # RUN RAT
    # ========================================================

    log "Running CAT_pack reads (RAT prok)"

    set +e

    "$CAT_CMD" reads \
        --mode cr \
        -c "$CONTIGS" \
        --c2c "$CAT_C2C" \
        --u2c "$DUMMY_U2C" \
        -1 "$R1" \
        -2 "$R2" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --tmpdir "$RAT_TMPDIR" \
        --read2classification \
        --force \
        --nproc "$THREADS" \
        2>&1 | tee "$LOGFILE"

    EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "RAT exit code: $EXIT_CODE"

    # ========================================================
    # HANDLE FAILURE
    # ========================================================

    if [[ "$EXIT_CODE" -ne 0 ]]; then

        log "ERROR: RAT prok failed for $SAMPLE"
        log "Last lines of RAT log:"

        tail -n 50 "$LOGFILE" || true

        rm -f "$DUMMY_U2C"

        return 1
    fi

    # ========================================================
    # VALIDATE OUTPUT
    # ========================================================

    [[ -s "$READ_FILE" ]] || {
        log "ERROR: Missing RAT read2classification output"
        rm -f "$DUMMY_U2C"
        return 1
    }

    [[ -s "$COMPLETE_FILE" ]] || {
        log "ERROR: Missing RAT complete abundance output"
        rm -f "$DUMMY_U2C"
        return 1
    }

    [[ -s "$CONTIG_FILE" ]] || {
        log "ERROR: Missing RAT contig abundance output"
        rm -f "$DUMMY_U2C"
        return 1
    }

    # ========================================================
    # SUMMARY
    # ========================================================

    log "RAT prok outputs generated:"

    ls -lh "$OUTDIR"/*.RAT.* 2>/dev/null || true

    # ========================================================
    # CLEANUP
    # ========================================================

    rm -f "$DUMMY_U2C"

    if [[ "${DEBUG:-0}" != "1" ]]; then

        rm -rf "$RAT_TMPDIR"

    else

        log "DEBUG=1 -> keeping RAT temporary directory:"
        log "$RAT_TMPDIR"

    fi

    log "Finished RAT prok for $SAMPLE"

    return 0
}

# ============================================================
# SAMPLE SELECTION
# ============================================================

# Backward compatibility:
# TARGET_SAMPLE still works, but SAMPLE is preferred.
if [[ -z "${SAMPLE:-}" && -n "${TARGET_SAMPLE:-}" ]]; then
    SAMPLE="$TARGET_SAMPLE"
fi

if [[ -n "${SAMPLE:-}" ]]; then

    # --------------------------------------------------------
    # SINGLE-SAMPLE MODE
    # --------------------------------------------------------

    log "Running in single-sample mode"

    run_rat_prok "$SAMPLE"

else

    # --------------------------------------------------------
    # ARRAY MODE
    # --------------------------------------------------------

    log "Running in array mode"

    # Stable sample list based on Tiara directories.
    # Samples with empty prok.fa simply skip inside
    # run_rat_prok().
    mapfile -t SAMPLE_DIRS < <(
        find "$TIARA_DIR" \
            -mindepth 1 \
            -maxdepth 1 \
            -type d \
            | sort
    )

    N=${#SAMPLE_DIRS[@]}

    log "Total Tiara sample directories detected: $N"

    [[ "$N" -gt 0 ]] || {
        log "ERROR: No Tiara sample directories found"
        exit 1
    }

    # -------------------- ARRAY CHECK --------------------

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Run as array job or export SAMPLE=<sample>"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} out of range (1-$N)"
        exit 1
    fi

    # -------------------- SELECT SAMPLE --------------------

    SAMPLE_DIR="${SAMPLE_DIRS[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$SAMPLE_DIR")

    log "Array task: ${SLURM_ARRAY_TASK_ID} / $N"

    run_rat_prok "$SAMPLE"

fi

log "RAT prokaryotic read-level step completed"

;;


# ============================================================
### RAT SUMMARY (STRUCTURED TAXONOMIC RANKS — READ LEVEL)
# ============================================================
rat_reads_summary_structured)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting RAT structured summary step"

RAT_DIR="$OUTPUT_DIR/rat_reads"
SUMMARY_DIR="$OUTPUT_DIR/rat_reads_summary"
LOG_DIR="$OUTPUT_DIR/logs/rat_reads_summary"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

NAMES="$CAT_TAX/names.dmp"
NODES="$CAT_TAX/nodes.dmp"

FORCE="${FORCE:-false}"

mkdir -p "$SUMMARY_DIR" "$LOG_DIR"

# ============================================================
# CHECK FILES
# ============================================================
[[ -s "$NAMES" ]] || {
    log "ERROR: names.dmp missing"
    exit 1
}

[[ -s "$NODES" ]] || {
    log "ERROR: nodes.dmp missing"
    exit 1
}

# ============================================================
# BUILD LOOKUP TABLES
# ============================================================

if [[ ! -s "$SUMMARY_DIR/taxid2name.tsv" ]]; then

    log "Building taxid2name.tsv"

    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")

        if ($4 == "scientific name") {
            print $1 "\t" $2
        }
    }
    ' "$NAMES" \
    > "$SUMMARY_DIR/taxid2name.tsv"
fi

if [[ ! -s "$SUMMARY_DIR/taxid2rank.tsv" ]]; then

    log "Building taxid2rank.tsv"

    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")
        print $1 "\t" $3
    }
    ' "$NODES" \
    > "$SUMMARY_DIR/taxid2rank.tsv"
fi

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================

for SAMPLE_DIR in "$RAT_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")

    OUT="$SUMMARY_DIR/${SAMPLE}.taxonomy_structured.tsv"
    TMP="$OUT.tmp"

    log "----------------------------------------"
    log "Processing $SAMPLE"

    # ========================================================
    # SKIP IF EXISTS
    # ========================================================
    if [[ -s "$OUT" && "$FORCE" != "true" ]]; then
        log "Skipping $SAMPLE (already processed)"
        continue
    fi

    # ========================================================
    # FIND RAT OUTPUT
    # ========================================================
    IN_FILES=("$SAMPLE_DIR"/*.RAT.contig.abundance.txt)

    if [[ ${#IN_FILES[@]} -eq 0 ]]; then
        log "Skipping $SAMPLE (no RAT output found)"
        continue
    fi

    IN="${IN_FILES[0]}"

    [[ -s "$IN" ]] || {
        log "Skipping $SAMPLE (input empty)"
        continue
    }

    log "Using input: $(basename "$IN")"

    # ========================================================
    # PROCESS RAT OUTPUT
    # ========================================================
    awk -F'\t' \
        -v names="$SUMMARY_DIR/taxid2name.tsv" \
        -v ranks="$SUMMARY_DIR/taxid2rank.tsv" '
    BEGIN {

        while ((getline < names) > 0) {
            name[$1] = $2
        }
        close(names)

        while ((getline < ranks) > 0) {
            rank[$1] = $2
        }
        close(ranks)

        OFS="\t"

        print "contig","reads","fraction","coverage","superkingdom","kingdom","phylum","class","order","family","genus","species"
    }

    /^#/ { next }

    $1 == "unmapped" { next }

    {

        contig = $1
        reads = $2
        fraction = $3
        coverage = $5
        lineage = $6

        if (lineage == "" || lineage == "no taxid assigned") {

            print contig,reads,fraction,coverage,"NA","NA","NA","NA","NA","NA","NA","NA"
            next
        }

        split(lineage, taxids, ";")

        sk="NA"
        k="NA"
        p="NA"
        c="NA"
        o="NA"
        f="NA"
        g="NA"
        s="NA"

        for (i in taxids) {

            t = taxids[i]

            r = rank[t]
            n = name[t]

            if (r == "superkingdom") sk = n
            else if (r == "kingdom") k = n
            else if (r == "phylum") p = n
            else if (r == "class") c = n
            else if (r == "order") o = n
            else if (r == "family") f = n
            else if (r == "genus") g = n
            else if (r == "species") s = n
        }

        print contig,reads,fraction,coverage,sk,k,p,c,o,f,g,s
    }
    ' "$IN" \
    > "$TMP"

    # ========================================================
    # SAFE WRITE
    # ========================================================
    if [[ -s "$TMP" ]]; then

        mv "$TMP" "$OUT"

        log "Created: $OUT"

    else

        log "WARNING: Empty output for $SAMPLE"

        rm -f "$TMP"
    fi

done

log "Finished RAT structured summary"
;;

# ============================================================
### RAT SUMMARY (STRUCTURED TAXONOMIC RANKS — READ LEVEL, PROK)
# ============================================================
rat_reads_summary_structured_prok)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting RAT structured summary step (prok)"

RAT_DIR="$OUTPUT_DIR/rat_reads_prok"

SUMMARY_DIR="$OUTPUT_DIR/rat_reads_summary_prok"

LOG_DIR="$OUTPUT_DIR/logs/rat_reads_summary_prok"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

NAMES="$CAT_TAX/names.dmp"

NODES="$CAT_TAX/nodes.dmp"

FORCE="${FORCE:-false}"

mkdir -p "$SUMMARY_DIR" "$LOG_DIR"

# ============================================================
# CHECK FILES
# ============================================================
[[ -s "$NAMES" ]] || {
    log "ERROR: names.dmp missing"
    exit 1
}

[[ -s "$NODES" ]] || {
    log "ERROR: nodes.dmp missing"
    exit 1
}

# ============================================================
# BUILD LOOKUP TABLES
# ============================================================

if [[ ! -s "$SUMMARY_DIR/taxid2name.tsv" ]]; then

    log "Building taxid2name.tsv"

    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")

        if ($4 == "scientific name") {
            print $1 "\t" $2
        }
    }
    ' "$NAMES" \
    > "$SUMMARY_DIR/taxid2name.tsv"
fi

if [[ ! -s "$SUMMARY_DIR/taxid2rank.tsv" ]]; then

    log "Building taxid2rank.tsv"

    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")
        print $1 "\t" $3
    }
    ' "$NODES" \
    > "$SUMMARY_DIR/taxid2rank.tsv"
fi

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================

for SAMPLE_DIR in "$RAT_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")

    OUT="$SUMMARY_DIR/${SAMPLE}.taxonomy_structured.tsv"

    TMP="$OUT.tmp"

    log "----------------------------------------"
    log "Processing $SAMPLE"

    # ========================================================
    # SKIP IF EXISTS
    # ========================================================
    if [[ -s "$OUT" && "$FORCE" != "true" ]]; then

        log "Skipping $SAMPLE (already processed)"

        continue
    fi

    # ========================================================
    # FIND RAT OUTPUT
    # ========================================================
    IN_FILES=("$SAMPLE_DIR"/*.RAT.read2classification.txt)

    if [[ ${#IN_FILES[@]} -eq 0 ]]; then

        log "Skipping $SAMPLE (no RAT output found)"

        continue
    fi

    IN="${IN_FILES[0]}"

    [[ -s "$IN" ]] || {

        log "Skipping $SAMPLE (input empty)"

        continue
    }

    log "Using input: $(basename "$IN")"

    # ========================================================
    # PROCESS RAT OUTPUT
    # ========================================================
    awk -F'\t' \
        -v names="$SUMMARY_DIR/taxid2name.tsv" \
        -v ranks="$SUMMARY_DIR/taxid2rank.tsv" '
    BEGIN {

        while ((getline < names) > 0) {
            name[$1] = $2
        }
        close(names)

        while ((getline < ranks) > 0) {
            rank[$1] = $2
        }
        close(ranks)

        OFS="\t"

        print "read","classification","superkingdom","kingdom","phylum","class","order","family","genus","species"
    }

    /^#/ { next }

    NR==1 { next }

    {

        read = $1
        lineage = $2

        if (lineage == "" || lineage == "no taxid assigned") {

            print read,"unclassified","NA","NA","NA","NA","NA","NA","NA","NA"

            next
        }

        split(lineage, taxids, ";")

        sk="NA"
        k="NA"
        p="NA"
        c="NA"
        o="NA"
        f="NA"
        g="NA"
        s="NA"

        for (i in taxids) {

            t = taxids[i]

            r = rank[t]
            n = name[t]

            if (r == "superkingdom") sk = n
            else if (r == "kingdom") k = n
            else if (r == "phylum") p = n
            else if (r == "class") c = n
            else if (r == "order") o = n
            else if (r == "family") f = n
            else if (r == "genus") g = n
            else if (r == "species") s = n
        }

        print read,"classified",sk,k,p,c,o,f,g,s
    }
    ' "$IN" \
    > "$TMP"

    # ========================================================
    # SAFE WRITE
    # ========================================================
    if [[ -s "$TMP" ]]; then

        mv "$TMP" "$OUT"

        log "Created: $OUT"

    else

        log "WARNING: Empty output for $SAMPLE"

        rm -f "$TMP"
    fi

done

log "Finished RAT structured summary (prok)"
;;

# ============================================================
### CAT SUMMARY (STRUCTURED TAXONOMIC RANKS — CONTIG LEVEL)
# ============================================================
cat_contigs_summary_structured)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting CAT structured summary step"

CAT_DIR="$OUTPUT_DIR/cat_contigs"
SUMMARY_DIR="$OUTPUT_DIR/cat_contigs_summary"
LOG_DIR="$OUTPUT_DIR/logs/cat_contigs_summary"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

NAMES="$CAT_TAX/names.dmp"
NODES="$CAT_TAX/nodes.dmp"

# Optional: allow overwrite
FORCE="${FORCE:-false}"

mkdir -p "$SUMMARY_DIR" "$LOG_DIR"

# -------------------- CHECK --------------------
[[ -s "$NAMES" ]] || { log "ERROR: names.dmp missing"; exit 1; }
[[ -s "$NODES" ]] || { log "ERROR: nodes.dmp missing"; exit 1; }

# ============================================================
# BUILD LOOKUP TABLES (REUSE IF EXISTS)
# ============================================================

if [[ ! -s "$SUMMARY_DIR/taxid2name.tsv" ]]; then
    log "Building taxid2name.tsv"
    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")
        if ($4 == "scientific name") {
            print $1 "\t" $2
        }
    }
    ' "$NAMES" > "$SUMMARY_DIR/taxid2name.tsv"
fi

if [[ ! -s "$SUMMARY_DIR/taxid2rank.tsv" ]]; then
    log "Building taxid2rank.tsv"
    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")
        print $1 "\t" $3
    }
    ' "$NODES" > "$SUMMARY_DIR/taxid2rank.tsv"
fi

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================

for SAMPLE_DIR in "$CAT_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")
    OUT="$SUMMARY_DIR/${SAMPLE}.taxonomy_structured.tsv"
    TMP="$OUT.tmp"

    log "----------------------------------------"
    log "Processing $SAMPLE"

    # -------------------- SKIP IF EXISTS --------------------
    if [[ -s "$OUT" && "$FORCE" != "true" ]]; then
        log "Skipping $SAMPLE (already processed)"
        continue
    fi

    # -------------------- FIND CAT OUTPUT (SAFE) --------------------
    IN_FILES=("$SAMPLE_DIR"/*contig*classification.txt)

    if [[ ${#IN_FILES[@]} -eq 0 ]]; then
        log "Skipping $SAMPLE (no CAT output found)"
        continue
    fi

    IN="${IN_FILES[0]}"

    [[ -s "$IN" ]] || {
        log "Skipping $SAMPLE (input empty)"
        continue
    }

    log "Using input: $(basename "$IN")"

    # -------------------- PROCESS --------------------
    awk -F'\t' \
        -v names="$SUMMARY_DIR/taxid2name.tsv" \
        -v ranks="$SUMMARY_DIR/taxid2rank.tsv" '
    BEGIN {
        # load names
        while ((getline < names) > 0) {
            name[$1] = $2
        }
        close(names)

        # load ranks
        while ((getline < ranks) > 0) {
            rank[$1] = $2
        }
        close(ranks)

        OFS="\t"

        print "contig","superkingdom","kingdom","phylum","class","order","family","genus","species"
    }

    NR>1 {

        contig = $1
        lineage = $4

        # handle unclassified
        if (lineage == "" || lineage == "no taxid assigned") {
            print contig,"NA","NA","NA","NA","NA","NA","NA","NA"
            next
        }

        split(lineage, taxids, ";")

        sk="NA"; k="NA"; p="NA"; c="NA"; o="NA"; f="NA"; g="NA"; s="NA"

        for (i in taxids) {
            t = taxids[i]

            r = rank[t]
            n = name[t]

            if (r == "superkingdom") sk = n
            else if (r == "kingdom") k = n
            else if (r == "phylum") p = n
            else if (r == "class") c = n
            else if (r == "order") o = n
            else if (r == "family") f = n
            else if (r == "genus") g = n
            else if (r == "species") s = n
        }

        print contig, sk, k, p, c, o, f, g, s
    }
    ' "$IN" > "$TMP"

    # -------------------- SAFE WRITE --------------------
    if [[ -s "$TMP" ]]; then
        mv "$TMP" "$OUT"
        log "Created: $OUT"
    else
        log "WARNING: Empty output for $SAMPLE"
        rm -f "$TMP"
    fi

done

log "Finished CAT structured summary"
;;

# ============================================================
### CAT SUMMARY (STRUCTURED TAXONOMIC RANKS — MITO CONTIGS)
# ============================================================
cat_mito_summary_structured)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting CAT MITO structured summary step"

CAT_DIR="$OUTPUT_DIR/cat_mito"
SUMMARY_DIR="$OUTPUT_DIR/cat_mito_summary"
LOG_DIR="$OUTPUT_DIR/logs/cat_mito_summary"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

NAMES="$CAT_TAX/names.dmp"
NODES="$CAT_TAX/nodes.dmp"

mkdir -p "$SUMMARY_DIR" "$LOG_DIR"

# -------------------- CHECK --------------------
[[ -s "$NAMES" ]] || { log "ERROR: names.dmp missing"; exit 1; }
[[ -s "$NODES" ]] || { log "ERROR: nodes.dmp missing"; exit 1; }

# ============================================================
# BUILD LOOKUP TABLES (REUSE IF EXISTS)
# ============================================================

if [[ ! -s "$SUMMARY_DIR/taxid2name.tsv" ]]; then
    log "Building taxid2name.tsv"
    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")
        if ($4 == "scientific name") {
            print $1 "\t" $2
        }
    }' "$NAMES" > "$SUMMARY_DIR/taxid2name.tsv"
fi

if [[ ! -s "$SUMMARY_DIR/taxid2rank.tsv" ]]; then
    log "Building taxid2rank.tsv"
    awk -F"\t\\|\t" '
    {
        gsub("\t\\|","")
        print $1 "\t" $3
    }' "$NODES" > "$SUMMARY_DIR/taxid2rank.tsv"
fi

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================

for SAMPLE_DIR in "$CAT_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")
    OUT="$SUMMARY_DIR/${SAMPLE}.mito.taxonomy_structured.tsv"

    log "----------------------------------------"
    log "Processing $SAMPLE"

    # -------------------- FIND CAT OUTPUT --------------------
    IN=$(ls "$SAMPLE_DIR"/*mito*classification.txt 2>/dev/null | head -n 1)

    [[ -s "$IN" ]] || {
        log "Skipping $SAMPLE (no mito CAT output found)"
        continue
    }

    log "Using input: $(basename "$IN")"

    # -------------------- PROCESS --------------------
    awk -F'\t' \
        -v names="$SUMMARY_DIR/taxid2name.tsv" \
        -v ranks="$SUMMARY_DIR/taxid2rank.tsv" '
    BEGIN {
        while ((getline < names) > 0) name[$1] = $2
        close(names)

        while ((getline < ranks) > 0) rank[$1] = $2
        close(ranks)

        OFS="\t"
        print "contig","superkingdom","kingdom","phylum","class","order","family","genus","species"
    }

    NR>1 {

        contig = $1
        lineage = $4

        if (lineage == "" || lineage == "no taxid assigned") {
            print contig,"NA","NA","NA","NA","NA","NA","NA","NA"
            next
        }

        split(lineage, taxids, ";")

        sk="NA"; k="NA"; p="NA"; c="NA"; o="NA"; f="NA"; g="NA"; s="NA"

        for (i in taxids) {
            t = taxids[i]
            r = rank[t]
            n = name[t]

            if (r == "superkingdom") sk = n
            else if (r == "kingdom") k = n
            else if (r == "phylum") p = n
            else if (r == "class") c = n
            else if (r == "order") o = n
            else if (r == "family") f = n
            else if (r == "genus") g = n
            else if (r == "species") s = n
        }

        print contig, sk, k, p, c, o, f, g, s
    }
    ' "$IN" > "$OUT"

    log "Created: $OUT"

done

log "Finished CAT MITO structured summary"
;;

# ============================================================
### BAT SUMMARY (STRUCTURED TAXONOMIC RANKS)
# ============================================================
bat_bins_summary_structured)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting BAT structured summary step"

BAT_DIR="$OUTPUT_DIR/bat_bins"
SUMMARY_DIR="$OUTPUT_DIR/bat_bins_summary"
LOG_DIR="$OUTPUT_DIR/logs/bat_bins_summary"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

NAMES="$CAT_TAX/names.dmp"
NODES="$CAT_TAX/nodes.dmp"

mkdir -p "$SUMMARY_DIR" "$LOG_DIR"

# -------------------- CHECK --------------------
[[ -s "$NAMES" ]] || { log "ERROR: names.dmp missing"; exit 1; }
[[ -s "$NODES" ]] || { log "ERROR: nodes.dmp missing"; exit 1; }

# ============================================================
# BUILD LOOKUP TABLES
# ============================================================

log "Preparing taxonomy lookup tables"

awk -F"\t\\|\t" '
{
    gsub("\t\\|","")
    if ($4 == "scientific name") {
        print $1 "\t" $2
    }
}
' "$NAMES" > "$SUMMARY_DIR/taxid2name.tsv"

awk -F"\t\\|\t" '
{
    gsub("\t\\|","")
    print $1 "\t" $3
}
' "$NODES" > "$SUMMARY_DIR/taxid2rank.tsv"

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================

for SAMPLE_DIR in "$BAT_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")
    IN="$SAMPLE_DIR/${SAMPLE}.bin2classification.txt"
    OUT="$SUMMARY_DIR/${SAMPLE}.taxonomy_structured.tsv"

    log "Processing $SAMPLE"

    [[ -s "$IN" ]] || {
        log "Skipping $SAMPLE (no BAT output)"
        continue
    }

    awk -F'\t' \
        -v names="$SUMMARY_DIR/taxid2name.tsv" \
        -v ranks="$SUMMARY_DIR/taxid2rank.tsv" '
    BEGIN {
        # load names
        while ((getline < names) > 0) {
            name[$1] = $2
        }
        close(names)

        # load ranks
        while ((getline < ranks) > 0) {
            rank[$1] = $2
        }
        close(ranks)

        OFS="\t"

        print "bin","superkingdom","kingdom","phylum","class","order","family","genus","species"
    }

    NR>1 {

        bin = $1
        lineage = $4

        split(lineage, taxids, ";")

        # initialize
        sk="NA"; k="NA"; p="NA"; c="NA"; o="NA"; f="NA"; g="NA"; s="NA"

        for (i in taxids) {
            t = taxids[i]

            r = rank[t]
            n = name[t]

            if (r == "superkingdom") sk = n
            else if (r == "kingdom") k = n
            else if (r == "phylum") p = n
            else if (r == "class") c = n
            else if (r == "order") o = n
            else if (r == "family") f = n
            else if (r == "genus") g = n
            else if (r == "species") s = n
        }

        print bin, sk, k, p, c, o, f, g, s
    }
    ' "$IN" > "$OUT"

    log "Created: $OUT"

done

log "Finished structured BAT summary"
;;

# ============================================================
### BAT SUMMARY (STRUCTURED TAXONOMIC RANKS — PROK MAGS)
# ============================================================
bat_prok_summary_structured)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting BAT prok structured summary step"

BAT_DIR="$OUTPUT_DIR/bat_prok"

SUMMARY_DIR="$OUTPUT_DIR/bat_prok_summary"

LOG_DIR="$OUTPUT_DIR/logs/bat_prok_summary"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

NAMES="$CAT_TAX/names.dmp"

NODES="$CAT_TAX/nodes.dmp"

mkdir -p "$SUMMARY_DIR" "$LOG_DIR"

# ============================================================
# CHECK TAXONOMY FILES
# ============================================================

[[ -s "$NAMES" ]] || {
    log "ERROR: names.dmp missing"
    exit 1
}

[[ -s "$NODES" ]] || {
    log "ERROR: nodes.dmp missing"
    exit 1
}

# ============================================================
# BUILD LOOKUP TABLES
# ============================================================

log "Preparing taxonomy lookup tables"

awk -F"\t\\|\t" '
{
    gsub("\t\\|","")

    if ($4 == "scientific name") {
        print $1 "\t" $2
    }
}
' "$NAMES" > "$SUMMARY_DIR/taxid2name.tsv"

awk -F"\t\\|\t" '
{
    gsub("\t\\|","")
    print $1 "\t" $3
}
' "$NODES" > "$SUMMARY_DIR/taxid2rank.tsv"

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================

N_DONE=0

for SAMPLE_DIR in "$BAT_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")

    IN="$SAMPLE_DIR/${SAMPLE}.bin2classification.txt"

    OUT="$SUMMARY_DIR/${SAMPLE}.taxonomy_structured.tsv"

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"

    # --------------------------------------------------------
    # CHECK INPUT
    # --------------------------------------------------------
    [[ -s "$IN" ]] || {

        log "Skipping $SAMPLE (missing BAT classification)"

        continue
    }

    # --------------------------------------------------------
    # PARSE BAT OUTPUT
    # --------------------------------------------------------
    awk -F'\t' \
        -v names="$SUMMARY_DIR/taxid2name.tsv" \
        -v ranks="$SUMMARY_DIR/taxid2rank.tsv" '
    BEGIN {

        # ----------------------------------------------------
        # LOAD SCIENTIFIC NAMES
        # ----------------------------------------------------
        while ((getline < names) > 0) {
            name[$1] = $2
        }
        close(names)

        # ----------------------------------------------------
        # LOAD TAXONOMIC RANKS
        # ----------------------------------------------------
        while ((getline < ranks) > 0) {
            rank[$1] = $2
        }
        close(ranks)

        OFS="\t"

        print "bin","superkingdom","kingdom","phylum","class","order","family","genus","species"
    }

    NR > 1 {

        bin = $1

        lineage = $4

        # remove spaces
        gsub(/ /, "", lineage)

        split(lineage, taxids, ";")

        # initialize
        sk="NA"
        k="NA"
        p="NA"
        c="NA"
        o="NA"
        f="NA"
        g="NA"
        s="NA"

        # ----------------------------------------------------
        # ASSIGN TAXONOMIC RANKS
        # ----------------------------------------------------
        for (i in taxids) {

            t = taxids[i]

            r = rank[t]

            n = name[t]

            if (r == "superkingdom") sk = n

            else if (r == "kingdom") k = n

            else if (r == "phylum") p = n

            else if (r == "class") c = n

            else if (r == "order") o = n

            else if (r == "family") f = n

            else if (r == "genus") g = n

            else if (r == "species") s = n
        }

        print bin, sk, k, p, c, o, f, g, s
    }
    ' "$IN" > "$OUT"

    # --------------------------------------------------------
    # VALIDATE OUTPUT
    # --------------------------------------------------------
    if [[ ! -s "$OUT" ]]; then

        log "ERROR: Failed creating summary for $SAMPLE"

        continue
    fi

    N_LINES=$(wc -l < "$OUT")

    log "Created: $OUT"
    log "Entries written: $((N_LINES - 1))"

    ((++N_DONE))

done

# ============================================================
# GLOBAL SUMMARY
# ============================================================

log "----------------------------------------"
log "Samples processed successfully: $N_DONE"

# ============================================================
# MERGE ALL SAMPLES
# ============================================================

MERGED="$SUMMARY_DIR/all_samples.taxonomy_structured.tsv"
MERGED_TMP=$(mktemp "$SUMMARY_DIR/.all_samples.taxonomy_structured.XXXXXX")

log "Creating merged summary table"

FIRST=1

for FILE in "$SUMMARY_DIR"/*.taxonomy_structured.tsv; do

    [[ "$FILE" == "$MERGED" ]] && continue
    [[ -s "$FILE" ]] || continue

    SAMPLE=$(basename "$FILE" .taxonomy_structured.tsv)

    if [[ "$FIRST" -eq 1 ]]; then

        awk -v sample="$SAMPLE" '
        BEGIN { OFS="\t" }

        NR==1 {
            print "sample",$0
            next
        }

        {
            print sample,$0
        }
        ' "$FILE" > "$MERGED_TMP"

        FIRST=0

    else

        awk -v sample="$SAMPLE" '
        BEGIN { OFS="\t" }

        NR>1 {
            print sample,$0
        }
        ' "$FILE" >> "$MERGED_TMP"
    fi

done

if [[ "$FIRST" -eq 1 ]]; then
    printf 'sample\tbin\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n' > "$MERGED_TMP"
fi
mv "$MERGED_TMP" "$MERGED"

if [[ -s "$MERGED" ]]; then

    TOTAL=$(($(wc -l < "$MERGED") - 1))

    log "Merged summary created: $MERGED"
    log "Total classified bins: $TOTAL"

else

    log "WARNING: merged summary not created"
fi

log "Finished BAT prok structured summary step"

;;

# ============================================================
### CONTIG TAXON ABUNDANCE (STRICT LEVEL, COVERAGE-WEIGHTED)
# ============================================================
cat_contigs_taxon_abundance)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

LEVEL="${LEVEL:-superkingdom}"

log "Running taxon abundance at level: $LEVEL"

OUT="$OUTPUT_DIR/cat_contigs_summary/taxon_abundance_${LEVEL}.tsv"

# -------------------- HEADER --------------------
echo -e "sample\ttaxon\tabundance" > "$OUT"

# ============================================================
# LOOP SAMPLES
# ============================================================
for FILE in "$OUTPUT_DIR/cat_contigs_summary"/*.taxonomy_structured.tsv; do

    SAMPLE=$(basename "$FILE" .taxonomy_structured.tsv)
    COV_FILE="$OUTPUT_DIR/coverage_coverm/${SAMPLE}.coverage.tsv"

    [[ -s "$COV_FILE" ]] || {
        log "Skipping $SAMPLE (no coverage file)"
        continue
    }

    log "Processing $SAMPLE"

    awk -F'\t' -v sample="$SAMPLE" -v level="$LEVEL" '

    BEGIN {
        OFS="\t"
    }

    # ========================================================
    # LOAD COVERAGE (coverage_weighted)
    # ========================================================
    FNR==NR {
        cov[$1] = $5
        next
    }

    NR==1 { next }

    {
        contig=$1

        sk=$2; k=$3; p=$4; c=$5; o=$6; f=$7; g=$8; s=$9

        # ====================================================
        # SELECT TAXONOMIC LEVEL
        # ====================================================
        if (level=="superkingdom") tax=sk
        else if (level=="kingdom") tax=k
        else if (level=="phylum") tax=p
        else if (level=="class") tax=c
        else if (level=="order") tax=o
        else if (level=="family") tax=f
        else if (level=="genus") tax=g
        else if (level=="species") tax=s
        else tax="NA"

        # ====================================================
        # STRICT FILTERING
        # ====================================================
        if (tax=="" || tax=="NA") next

        # remove literal artifacts (INCLUDING superkingdom)
        if (tax ~ /^(superkingdom|kingdom|phylum|class|order|family|genus|species)$/) next

        # ====================================================
        # VALIDATE KINGDOM LEVEL (very important)
        # ====================================================
        if (level=="kingdom") {
            if (!(tax=="Bacteria" || tax=="Archaea" || tax=="Eukaryota" || tax=="Viruses")) next
        }

        # ====================================================
        # COVERAGE
        # ====================================================
        val = cov[contig]

        if (val=="" || val==0) next

        sum[tax] += val
    }

    END {
        for (t in sum) {
            print sample, t, sum[t]
        }
    }

    ' "$COV_FILE" "$FILE" >> "$OUT"

done

log "Created: $OUT"

;;


# ============================================================
# DEFAULT
# ============================================================
*)
log "Unknown step: $STEP"
exit 1
;;

esac

log "Completed step: $STEP"
done

log "Pipeline finished successfully"
