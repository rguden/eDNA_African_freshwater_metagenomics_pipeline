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
    STEPS=(read_counts)
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
mapfile -t R1_FILES < <(ls "$READS_DIR"/*_R1*.fastq.gz | sort)

N=${#R1_FILES[@]}

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
mapfile -t R1_FILES < <(ls "$DEDUP_DIR"/*_dedup_R1.fastq.gz | sort)
N=${#R1_FILES[@]}

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
# PER-SAMPLE ASSEMBLY (ARRAY-AWARE)
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
log "TRIM_DIR: $TRIM_DIR"

# -------------------- FILE LIST --------------------
mapfile -t R1_FILES < <(ls "$TRIM_DIR"/*_trimmed_R1.fastq.gz | sort)
N=${#R1_FILES[@]}

log "Total samples detected: $N"

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
R1="${R1_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
SAMPLE=$(basename "$R1" | sed -E 's/_trimmed_R1.*//')

OUTDIR="$ASM_DIR/$SAMPLE"
CONTIGS="$OUTDIR/final.contigs.fa"
FINAL="$OUTDIR/${SAMPLE}.contigs.fa"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# -------------------- INPUT CHECK --------------------
[[ -s "$R1" && -s "$R2" ]] || {
    log "ERROR: Missing reads:"
    log "R1: $R1"
    log "R2: $R2"
    exit 1
}

# -------------------- SKIP IF DONE --------------------
if [[ -s "$FINAL" && "${FORCE:-0}" != "1" ]]; then
    log "Assembly exists → skipping"
    exit 0
fi

# -------------------- CLEAN OLD RUN --------------------
if [[ -d "$OUTDIR" ]]; then
    log "Removing existing directory: $OUTDIR"
    rm -rf "$OUTDIR"
fi

# DO NOT create OUTDIR — MEGAHIT will do it

# -------------------- RUN MEGAHIT --------------------
log "Running MEGAHIT"

megahit \
    -1 "$R1" \
    -2 "$R2" \
    -t "$THREADS" \
    --min-contig-len "$MINLEN" \
    --mem-flag 0 \
    -o "$OUTDIR" \
    &>> "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
[[ -s "$CONTIGS" ]] || {
    log "ERROR: Assembly failed"
    exit 1
}

# -------------------- RENAME --------------------
TMP="$OUTDIR/${SAMPLE}.renamed.fa"

awk -v s="$SAMPLE" '
/^>/ {gsub(/^>/,"",$0); print ">" s "_" $0; next}
{print}
' "$CONTIGS" > "$TMP"

mv "$TMP" "$FINAL"

# Optional: keep original but link it
ln -sfn "$(basename "$FINAL")" "$CONTIGS"

log "Finished: $FINAL"
;;

# ============================================================
# PER-SAMPLE MAPPING (ARRAY-AWARE)
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

# -------------------- TMPDIR VALIDATION --------------------
if [[ -z "${TMPDIR:-}" || ! -d "$TMPDIR" ]]; then
    log "ERROR: TMPDIR not set or not a directory"
    exit 1
fi

# Test write permission
touch "$TMPDIR/.tmp_test_${SLURM_JOB_ID}" 2>/dev/null || {
    log "ERROR: TMPDIR not writable: $TMPDIR"
    exit 1
}
rm -f "$TMPDIR/.tmp_test_${SLURM_JOB_ID}"

# Check free space (require ~20 GB; adjust if needed)
FREE_GB=$(df -BG "$TMPDIR" | awk 'NR==2{gsub("G","",$4); print $4}')
MIN_GB=${TMP_MIN_GB:-20}
if [[ "$FREE_GB" -lt "$MIN_GB" ]]; then
    log "ERROR: TMPDIR has only ${FREE_GB}G free (< ${MIN_GB}G required)"
    exit 1
fi

# Create per-job working dir
WORKDIR="$TMPDIR/map_${SLURM_JOB_ID}"
mkdir -p "$WORKDIR"

log "Using TMPDIR: $TMPDIR"
log "WORKDIR: $WORKDIR (free ${FREE_GB}G)"

# -------------------- TOOL CHECK --------------------
command -v bbmap.sh >/dev/null || { log "ERROR: bbmap.sh not found"; exit 1; }
command -v samtools  >/dev/null || { log "ERROR: samtools not found"; exit 1; }

# -------------------- COLLECT ASSEMBLIES --------------------
mapfile -t ASM_FILES < <(
    for d in "$OUTPUT_DIR/assembly"/*; do
        [[ -d "$d" ]] || continue
        s=$(basename "$d")
        f="$d/${s}.contigs.fa"
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No valid assemblies found"; exit 1; }

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$ASM" | sed 's/.contigs.fa//')

R1="$OUTPUT_DIR/trimmed/${SAMPLE}_trimmed_R1.fastq.gz"
R2="$OUTPUT_DIR/trimmed/${SAMPLE}_trimmed_R2.fastq.gz"
BAM="$MAP_DIR/${SAMPLE}.bam"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "ASM: $ASM"
log "R1:  $R1"
log "R2:  $R2"

# -------------------- INPUT CHECK --------------------
[[ -s "$ASM" ]] || { log "ERROR: Assembly missing"; exit 1; }
[[ -s "$R1" && -s "$R2" ]] || {
    log "ERROR: Missing reads"
    log "R1 exists? $( [[ -s "$R1" ]] && echo YES || echo NO )"
    log "R2 exists? $( [[ -s "$R2" ]] && echo YES || echo NO )"
    exit 1
}

# -------------------- SKIP --------------------
if [[ -s "$BAM" && "${FORCE:-0}" != "1" ]]; then
    log "Mapping exists → skipping"
    exit 0
fi

# -------------------- TMP FILES --------------------
TMP_PREFIX="$WORKDIR/${SAMPLE}"
TMP_REF="${TMP_PREFIX}.fa"
TMP_R1="${TMP_PREFIX}_R1.fastq.gz"
TMP_R2="${TMP_PREFIX}_R2.fastq.gz"
TMP_BAM="${TMP_PREFIX}.bam"

log "Copying inputs to WORKDIR"
cp "$ASM" "$TMP_REF"
cp "$R1" "$TMP_R1"
cp "$R2" "$TMP_R2"

# -------------------- THREAD SPLIT --------------------
BBMAP_THREADS=$((THREADS - 4))
SORT_THREADS=4
(( BBMAP_THREADS < 1 )) && BBMAP_THREADS=1

# -------------------- MEMORY (dynamic) --------------------
MEM_MB="${SLURM_MEM_PER_NODE:-0}"
if [[ "$MEM_MB" -gt 0 ]]; then
    MEM_GB=$((MEM_MB / 1024))
    JAVA_MEM=$((MEM_GB * 80 / 100))
else
    JAVA_MEM=16
fi

log "BBMap threads: $BBMAP_THREADS"
log "Sort threads: $SORT_THREADS"
log "Java memory: ${JAVA_MEM}g"

# -------------------- RUN BBMAP --------------------
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
| samtools view -b - \
| samtools sort \
    -@ "$SORT_THREADS" \
    -T "$WORKDIR/${SAMPLE}" \
    -o "$TMP_BAM" \
    2>> "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
[[ -s "$TMP_BAM" ]] || { log "ERROR: BAM not created"; exit 1; }

# -------------------- COPY BACK --------------------
log "Copying BAM back to output"
cp "$TMP_BAM" "$BAM"
samtools index "$BAM"

# -------------------- CLEANUP --------------------
if [[ "${DEBUG:-0}" != "1" ]]; then
    rm -rf "$WORKDIR"
fi

log "Finished $SAMPLE"
;;

# ============================================================
# QUAST (ASSEMBLY QC — ARRAY-AWARE)
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

# -------------------- COLLECT ASSEMBLIES (STRICT, SAME AS MAPPING) --------------------
mapfile -t ASM_FILES < <(
    for d in "$ASM_DIR"/*; do
        [[ -d "$d" ]] || continue
        SAMPLE=$(basename "$d")
        ASM="$d/${SAMPLE}.contigs.fa"
        [[ -s "$ASM" ]] && echo "$ASM"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total valid assemblies: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No valid assemblies found"
    exit 1
}

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID > N )); then
    log "Skipping: Task ID ${SLURM_ARRAY_TASK_ID} > $N"
    exit 0
fi

# -------------------- SELECT SAMPLE --------------------
ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$ASM" | sed 's/.contigs.fa//')

OUTDIR="$QUAST_DIR/$SAMPLE"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "Assembly: $ASM"
log "Output: $OUTDIR"

# -------------------- INPUT CHECK --------------------
[[ -s "$ASM" ]] || {
    log "ERROR: Missing assembly"
    exit 1
}

# -------------------- SKIP --------------------
if [[ -s "$OUTDIR/report.txt" && "${FORCE:-0}" != "1" ]]; then
    log "QUAST exists → skipping"
    exit 0
fi

# -------------------- CLEAN --------------------
rm -rf "$OUTDIR"

# -------------------- RUN QUAST --------------------
log "Running QUAST"

quast.py \
    "$ASM" \
    -o "$OUTDIR" \
    -t "$THREADS" \
    --min-contig 1000 \
    &>> "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
if [[ ! -s "$OUTDIR/report.txt" ]]; then
    log "ERROR: QUAST failed for $SAMPLE"
    exit 1
fi

log "Finished QUAST: $OUTDIR/report.txt"

;;

# ============================================================
# QUAST SUMMARY (MERGE ALL SAMPLES INTO TSV)
# ============================================================
quast_summary)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

QUAST_DIR="$OUTPUT_DIR/quast"
META="/staging/leuven/stg_00100/Rodgee/caves_data/metadata/biological_metadata.tsv"
OUT="$OUTPUT_DIR/quast/quast_summary.tsv"

mkdir -p "$(dirname "$OUT")"

log "Creating QUAST summary table"

# -------------------- HEADER --------------------
echo -e "sample\tsample_code\tincubation\tn_contigs\ttotal_length\tN50\tGC" > "$OUT"

# -------------------- LOOP --------------------
for REPORT in "$QUAST_DIR"/*/report.txt; do

    [[ -s "$REPORT" ]] || continue

    SAMPLE=$(basename "$(dirname "$REPORT")")

    # -------------------- EXTRACT METADATA --------------------
    read INC CODE < <(
        awk -v s="$SAMPLE" '
        BEGIN {FS="\t"}
        $1==s {print $2, $3}
        ' "$META"
    )

    # -------------------- EXTRACT QUAST METRICS --------------------
    ncontigs=$(grep -P "^# contigs \(>= 0 bp\)" "$REPORT" | awk '{print $NF}')
    totlen=$(grep -P "^Total length \(>= 0 bp\)" "$REPORT" | awk '{print $NF}')
    n50=$(grep -P "^N50" "$REPORT" | awk '{print $NF}')
    gc=$(grep -P "^GC \(%\)" "$REPORT" | awk '{print $NF}')

    echo -e "${SAMPLE}\t${CODE}\t${INC}\t${ncontigs}\t${totlen}\t${n50}\t${gc}" >> "$OUT"

done

log "Created: $OUT"

;;

# ============================================================
# QUAST COMPARISON (SHORT vs LONG — FIXED)
# ============================================================
quast_compare_incubation)

set -euo pipefail

log() { echo "[$(date '+%F %T')] $*"; }

IN="$OUTPUT_DIR/quast/quast_summary.tsv"
OUT="$OUTPUT_DIR/quast/quast_comparison.tsv"

log "Creating QUAST paired comparison"

awk '
BEGIN {
    FS="\t"; OFS="\t"
}

NR==1 { next }

{
    # FIX: normalize sample_code (remove -L / -S)
    base = $2
    sub(/-[LS]$/, "", base)

    if ($3 == "long") {
        long[base] = $0
    } else if ($3 == "short") {
        short[base] = $0
    }
}

END {
    print "sample_code",
          "long_sample","short_sample",
          "long_contigs","short_contigs",
          "long_length","short_length",
          "long_N50","short_N50",
          "long_GC","short_GC",
          "delta_N50","delta_length"

    for (k in long) {
        if (k in short) {

            split(long[k], L, "\t")
            split(short[k], S, "\t")

            print k,
                  L[1], S[1],
                  L[4], S[4],
                  L[5], S[5],
                  L[6], S[6],
                  L[7], S[7],
                  (L[6] - S[6]),
                  (L[5] - S[5])
        }
    }
}
' "$IN" > "$OUT"

log "Created: $OUT"

;;

# ============================================================
# PER-SAMPLE COVERAGE (COVERM — ARRAY-AWARE)
# ============================================================
coverage_coverm)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

THREADS="${THREADS:-8}"

ASM_DIR="$OUTPUT_DIR/assembly"
MAP_DIR="$OUTPUT_DIR/mapping"
COV_DIR="$OUTPUT_DIR/coverage"
LOG_DIR="$OUTPUT_DIR/logs/coverage"

mkdir -p "$COV_DIR" "$LOG_DIR"

log "Starting CoverM step"
log "THREADS: $THREADS"

# -------------------- TOOL CHECK --------------------
command -v coverm >/dev/null || {
    log "ERROR: coverm not found (activate conda env)"
    exit 1
}

# -------------------- COLLECT ASSEMBLIES --------------------
mapfile -t ASM_FILES < <(
    for d in "$ASM_DIR"/*; do
        [[ -d "$d" ]] || continue
        SAMPLE=$(basename "$d")
        ASM="$d/${SAMPLE}.contigs.fa"
        [[ -s "$ASM" ]] && echo "$ASM"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No assemblies found"
    exit 1
}

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$ASM" | sed 's/.contigs.fa//')

BAM="$MAP_DIR/${SAMPLE}.bam"
RAW="$COV_DIR/${SAMPLE}.coverage.raw.tsv"
OUT="$COV_DIR/${SAMPLE}.coverage.tsv"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# -------------------- INPUT CHECK --------------------
[[ -s "$ASM" ]] || { log "ERROR: Missing assembly"; exit 1; }
[[ -s "$BAM" ]] || { log "ERROR: Missing BAM"; exit 1; }

# ============================================================
# STEP 1: RUN COVERM (ONLY IF RAW MISSING)
# ============================================================
if [[ -s "$RAW" && "${FORCE:-0}" != "1" ]]; then
    log "Raw coverage exists → skipping CoverM"
else
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

    [[ -s "$RAW" ]] || {
        log "ERROR: CoverM failed or empty output"
        exit 1
    }
fi

# ============================================================
# STEP 2: COMPUTE WEIGHTED COVERAGE (ROBUST PARSING)
# ============================================================
if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
    log "Final coverage exists → skipping weighting"
    exit 0
fi

log "Computing coverage_weighted"

# Debug header (important)
log "Header:"
head -n 1 "$RAW" | tee -a "$LOGFILE"

# -------------------- ROBUST PARSER --------------------
awk '
BEGIN {OFS="\t"}

NR==1 {
    print "contig", "mean_coverage", "covered_fraction", "length", "coverage_weighted"
    next
}

{
    contig = $1

    # collect ONLY numeric fields
    n = 0
    for (i=1; i<=NF; i++) {
        if ($i ~ /^[0-9.]+$/) {
            nums[++n] = $i
        }
    }

    if (n < 3) next

    cov  = nums[n-2]
    frac = nums[n-1]
    len  = nums[n]

    print contig, cov, frac, len, cov * len
}
' "$RAW" > "$OUT"

# -------------------- VALIDATE --------------------
[[ -s "$OUT" ]] || {
    log "ERROR: Output not created"
    exit 1
}

# -------------------- SUMMARY --------------------
log "Preview:"
head -n 5 "$OUT" | tee -a "$LOGFILE"

log "Finished coverage for $SAMPLE"

;;

# ============================================================
# COVERAGE SUMMARY PER SAMPLE
# ============================================================
coverage_summary_sample)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

COV_DIR="$OUTPUT_DIR/coverage"
OUT="$OUTPUT_DIR/coverage/coverage_per_sample.tsv"

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
# PER-SAMPLE COVERAGE (COVERM — ARRAY-AWARE)
# ============================================================
coverage_coverm)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

THREADS="${THREADS:-8}"

ASM_DIR="$OUTPUT_DIR/assembly"
MAP_DIR="$OUTPUT_DIR/mapping"
COV_DIR="$OUTPUT_DIR/coverage"
LOG_DIR="$OUTPUT_DIR/logs/coverage"

mkdir -p "$COV_DIR" "$LOG_DIR"

log "Starting CoverM step"
log "THREADS: $THREADS"

# -------------------- TOOL CHECK --------------------
command -v coverm >/dev/null || {
    log "ERROR: coverm not found (activate conda env)"
    exit 1
}

# -------------------- COLLECT ASSEMBLIES --------------------
mapfile -t ASM_FILES < <(
    for d in "$ASM_DIR"/*; do
        [[ -d "$d" ]] || continue
        SAMPLE=$(basename "$d")
        ASM="$d/${SAMPLE}.contigs.fa"
        [[ -s "$ASM" ]] && echo "$ASM"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No assemblies found"
    exit 1
}

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
ASM="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$ASM" | sed 's/.contigs.fa//')

BAM="$MAP_DIR/${SAMPLE}.bam"
RAW="$COV_DIR/${SAMPLE}.coverage.raw.tsv"
OUT="$COV_DIR/${SAMPLE}.coverage.tsv"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# -------------------- INPUT CHECK --------------------
[[ -s "$ASM" ]] || { log "ERROR: Missing assembly"; exit 1; }
[[ -s "$BAM" ]] || { log "ERROR: Missing BAM"; exit 1; }

# ============================================================
# STEP 1: RUN COVERM (ONLY IF RAW MISSING)
# ============================================================
if [[ -s "$RAW" && "${FORCE:-0}" != "1" ]]; then
    log "Raw coverage exists → skipping CoverM"
else
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

    [[ -s "$RAW" ]] || {
        log "ERROR: CoverM failed or empty output"
        exit 1
    }
fi

# ============================================================
# STEP 2: COMPUTE WEIGHTED COVERAGE (ROBUST PARSING)
# ============================================================
if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
    log "Final coverage exists → skipping weighting"
    exit 0
fi

log "Computing coverage_weighted"

# Debug header (important)
log "Header:"
head -n 1 "$RAW" | tee -a "$LOGFILE"

# -------------------- ROBUST PARSER --------------------
awk '
BEGIN {OFS="\t"}

NR==1 {
    print "contig", "mean_coverage", "covered_fraction", "length", "coverage_weighted"
    next
}

{
    contig = $1

    # collect ONLY numeric fields
    n = 0
    for (i=1; i<=NF; i++) {
        if ($i ~ /^[0-9.]+$/) {
            nums[++n] = $i
        }
    }

    if (n < 3) next

    cov  = nums[n-2]
    frac = nums[n-1]
    len  = nums[n]

    print contig, cov, frac, len, cov * len
}
' "$RAW" > "$OUT"

# -------------------- VALIDATE --------------------
[[ -s "$OUT" ]] || {
    log "ERROR: Output not created"
    exit 1
}

# -------------------- SUMMARY --------------------
log "Preview:"
head -n 5 "$OUT" | tee -a "$LOGFILE"

log "Finished coverage for $SAMPLE"

;;

# ============================================================
# PER-SAMPLE COVERAGE (SAMTOOLS -AA — ARRAY-AWARE)
# ============================================================
coverage_samtools)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

THREADS="${THREADS:-8}"

MAP_DIR="$OUTPUT_DIR/mapping"
COV_DIR="$OUTPUT_DIR/coverage_samtools"
LOG_DIR="$OUTPUT_DIR/logs/coverage_samtools"

mkdir -p "$COV_DIR" "$LOG_DIR"

log "Starting samtools coverage step"

# -------------------- COLLECT BAM FILES --------------------
mapfile -t BAM_FILES < <(
    for f in "$MAP_DIR"/*.bam; do
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#BAM_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No BAM files found"
    exit 1
}

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
BAM="${BAM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$BAM" .bam)

OUT="$COV_DIR/${SAMPLE}.samtools.tsv"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# -------------------- INPUT CHECK --------------------
[[ -s "$BAM" ]] || { log "ERROR: Missing BAM"; exit 1; }

# ============================================================
# STEP: SAMTOOLS DEPTH
# ============================================================
if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
    log "Output exists → skipping"
    exit 0
fi

log "Running samtools depth (-aa)"

samtools depth -aa "$BAM" | \
awk -v sample="$SAMPLE" '
BEGIN {
    sum = 0;
    total = 0;
    covered = 0;
}
{
    depth = $3;

    total++;

    if (depth > 0) {
        covered++;
        sum += depth;
    }
}
END {
    mean = (total > 0) ? sum / total : 0;
    breadth = (total > 0) ? covered / total : 0;

    print "sample\tmean_coverage\tbreadth";
    printf "%s\t%.6f\t%.6f\n", sample, mean, breadth;
}
' > "$OUT" 2>> "$LOGFILE"

# -------------------- VALIDATE --------------------
[[ -s "$OUT" ]] || {
    log "ERROR: samtools output not created"
    exit 1
}

# -------------------- PREVIEW --------------------
log "Preview:"
head -n 5 "$OUT" | tee -a "$LOGFILE"

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
### COVERAGE (JGI DEPTH — PER-SAMPLE, ARRAY-AWARE)
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

# -------------------- COLLECT SAMPLES --------------------
mapfile -t BAM_FILES < <(
    for f in "$MAP_DIR"/*.bam; do
        [[ -f "$f" ]] && echo "$f"
    done | sort
)

N=${#BAM_FILES[@]}
log "Total BAM files: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No BAM files found"; exit 1; }

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID > N )); then
    log "No sample for this index → exiting"
    exit 0
fi

# -------------------- SELECT SAMPLE --------------------
BAM="${BAM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$BAM" .bam)

DEPTH="$COV_DIR/${SAMPLE}.depth.txt"

log "----------------------------------------"
log "Processing: $SAMPLE"
log "BAM: $BAM"

# -------------------- SKIP --------------------
if [[ -s "$DEPTH" && "${FORCE:-0}" != "1" ]]; then
    log "Depth already exists → skipping"
    exit 0
fi

# -------------------- RUN JGI DEPTH --------------------
log "Running jgi_summarize_bam_contig_depths"

jgi_summarize_bam_contig_depths \
    --outputDepth "$DEPTH" \
    "$BAM" \
    &>> "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
[[ -s "$DEPTH" ]] || {
    log "ERROR: Depth calculation failed"
    exit 1
}

log "Finished depth for $SAMPLE"
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
# ============================================================
binning_metabat_prok)

set -euo pipefail
shopt -s nullglob

log "Starting MetaBAT2 binning (clean assembly)"

THREADS="${THREADS:-8}"
MIN_LEN="${MIN_LEN:-1500}"

ASM_DIR="$OUTPUT_DIR/assembly_clean"
COV_DIR="$OUTPUT_DIR/coverage_metabat"
COV_FILT_DIR="$OUTPUT_DIR/coverage_metabat_filtered"

BIN_DIR="$OUTPUT_DIR/binning_metabat_prok"
LOG_DIR="$OUTPUT_DIR/logs/binning_metabat_prok"

mkdir -p "$BIN_DIR" "$LOG_DIR"

# -------------------- COLLECT SAMPLES --------------------
mapfile -t ASM_FILES < <(
    for d in "$ASM_DIR"/*; do
        [[ -d "$d" ]] || continue
        s=$(basename "$d")
        f="$d/${s}.contigs.clean.fa"
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total assemblies: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No cleaned assemblies found"; exit 1; }

# ============================================================
# SINGLE-SAMPLE OVERRIDE
# ============================================================
if [[ -n "${TARGET_SAMPLE:-}" ]]; then

    log "TARGET_SAMPLE detected: $TARGET_SAMPLE"

    FOUND=0

    for f in "${ASM_FILES[@]}"; do

        s=$(basename "$f" .contigs.clean.fa)

        if [[ "$s" == "$TARGET_SAMPLE" ]]; then
            ASM_FILES=("$f")
            FOUND=1
            break
        fi
    done

    if [[ "$FOUND" -eq 0 ]]; then
        log "ERROR: TARGET_SAMPLE not found in assemblies"
        exit 1
    fi

    N=1
    SLURM_ARRAY_TASK_ID=1
fi

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID > N )); then
    log "No sample for this index → exiting"
    exit 0
fi

# -------------------- SELECT SAMPLE --------------------
CONTIGS="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$CONTIGS" .contigs.clean.fa)

DEPTH_FILT="$COV_FILT_DIR/${SAMPLE}.depth.txt"
DEPTH_RAW="$COV_DIR/${SAMPLE}.depth.txt"

OUTDIR="$BIN_DIR/$SAMPLE"

log "----------------------------------------"
log "Processing: $SAMPLE"

# -------------------- INPUT CHECK --------------------
[[ -s "$CONTIGS" ]] || { log "No contigs → skip"; exit 0; }

# Prefer filtered depth
if [[ -s "$DEPTH_FILT" ]]; then
    DEPTH="$DEPTH_FILT"
    log "Using FILTERED depth"
elif [[ -s "$DEPTH_RAW" ]]; then
    DEPTH="$DEPTH_RAW"
    log "WARNING: Using UNFILTERED depth"
else
    log "No depth file → skip"
    exit 0
fi

log "CONTIGS: $CONTIGS"
log "DEPTH:   $DEPTH"

# HANDLE INCOMPLETE OUTPUT

if [[ -d "$OUTDIR" ]] && ! compgen -G "$OUTDIR/bin*.fa" > /dev/null; then
    log "Incomplete run detected → cleaning"
    rm -rf "$OUTDIR"
fi

# -------------------- SKIP --------------------
if [[ -d "$OUTDIR" ]] && compgen -G "$OUTDIR/bin*.fa" > /dev/null \
   && [[ "${FORCE:-0}" != "1" ]]; then
    log "Binning already exists → skipping"
    exit 0
fi

# -------------------- CLEAN --------------------
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# FILTER CONTIGS BY LENGTH

FILTERED_CONTIGS="$OUTDIR/${SAMPLE}.min${MIN_LEN}.fa"

log "Filtering contigs (min length = $MIN_LEN)"

seqkit seq -m "$MIN_LEN" "$CONTIGS" > "$FILTERED_CONTIGS"

[[ -s "$FILTERED_CONTIGS" ]] || {
    log "ERROR: No contigs after filtering"
    exit 1
}

# CONSISTENCY CHECK

log "Checking contig–depth consistency"

grep -v '^contigName' "$DEPTH" | cut -f1 | sort > "$OUTDIR/depth.list"
grep "^>" "$FILTERED_CONTIGS" | sed 's/^>//' | sort > "$OUTDIR/contigs.list"

comm -12 "$OUTDIR/depth.list" "$OUTDIR/contigs.list" > "$OUTDIR/shared.list"

N_SHARED=$(wc -l < "$OUTDIR/shared.list")

if [[ "$N_SHARED" -lt 100 ]]; then
    log "WARNING: Low overlap between contigs and depth ($N_SHARED)"
else
    log "Matching contigs: $N_SHARED"
fi

# RUN METABAT2

log "Running MetaBAT2"

metabat2 \
    -i "$FILTERED_CONTIGS" \
    -a "$DEPTH" \
    -o "$OUTDIR/bin" \
    -t "$THREADS" \
    -m "$MIN_LEN" \
    &>> "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
N_BINS=$(ls "$OUTDIR"/bin*.fa 2>/dev/null | wc -l || true)

if [[ "$N_BINS" -eq 0 ]]; then
    log "WARNING: No bins produced"
else
    log "Bins produced: $N_BINS"
fi

log "Finished MetaBAT2 for $SAMPLE"
;;

# ============================================================
### COVERM (CONTIG COVERAGE FOR MAXBIN— ARRAY-AWARE, FULL ASSEMBLY)
# ============================================================
coverm_prok)

set -euo pipefail
shopt -s nullglob

log "Starting CoverM (full contig coverage)"

THREADS="${THREADS:-8}"

MAP_DIR="$OUTPUT_DIR/mapping"
COV_DIR="$OUTPUT_DIR/coverage_coverm"
LOG_DIR="$OUTPUT_DIR/logs/coverage_coverm"

mkdir -p "$COV_DIR" "$LOG_DIR"

# -------------------- COLLECT SAMPLES --------------------
mapfile -t BAM_FILES < <(
    for f in "$MAP_DIR"/*.bam; do
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#BAM_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No BAM files found"; exit 1; }

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "No sample → exiting"
    exit 0
fi

# -------------------- SELECT SAMPLE --------------------
BAM="${BAM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$BAM" .bam)

OUT="$COV_DIR/${SAMPLE}.coverm.tsv"
LOGFILE="$LOG_DIR/${SAMPLE}.log"

log "----------------------------------------"
log "Processing: $SAMPLE"
log "BAM: $BAM"

# -------------------- INPUT CHECK --------------------
[[ -s "$BAM" ]] || {
    log "Missing BAM → skipping"
    exit 0
}

# -------------------- SKIP --------------------
if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
    log "CoverM already done → skipping"
    exit 0
fi

# -------------------- RUN COVERM --------------------
log "Running CoverM"

coverm contig \
    --bam-files "$BAM" \
    --methods mean covered_fraction length \
    --threads "$THREADS" \
    --output-file "$OUT" \
    2> "$LOGFILE"

# -------------------- VALIDATE --------------------
[[ -s "$OUT" ]] || {
    log "ERROR: CoverM failed for $SAMPLE"
    exit 1
}

# -------------------- PREVIEW --------------------
log "Preview:"
head -n 5 "$OUT" | tee -a "$LOGFILE"

log "Finished CoverM for $SAMPLE"
;;

# ============================================================
### CLEAN COVERM (MAKE MAXBIN-COMPATIBLE)
# ============================================================
clean_coverm)

set -euo pipefail
shopt -s nullglob

log "Starting CoverM cleaning"

COV_DIR="$OUTPUT_DIR/coverage_coverm"
ASM_CLEAN_DIR="$OUTPUT_DIR/assembly_clean"
OUT_DIR="$OUTPUT_DIR/coverage_coverm_clean"
LOG_DIR="$OUTPUT_DIR/logs/coverage_coverm_clean"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# -------------------- COLLECT FILES --------------------
mapfile -t COV_FILES < <(
    for f in "$COV_DIR"/*.coverm.tsv; do
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#COV_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No CoverM files found"
    exit 1
}

# MODE SELECTION

if [[ -n "${TARGET_SAMPLE:-}" ]]; then
    log "Running in MANUAL mode ($TARGET_SAMPLE)"

    COV="$COV_DIR/${TARGET_SAMPLE}.coverm.tsv"
    SAMPLE="$TARGET_SAMPLE"
    process_sample=true

elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "Running in ARRAY mode"

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "ERROR: Task ID out of range"
        exit 1
    fi

    COV="${COV_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
    SAMPLE=$(basename "$COV" .coverm.tsv)
    process_sample=true

else
    log "Running in NON-ARRAY mode"
    process_sample=false
fi


# FUNCTION

run_clean_coverm() {

    local COV="$1"
    local SAMPLE="$2"

    CLEAN_CONTIGS="$ASM_CLEAN_DIR/$SAMPLE/${SAMPLE}.contigs.clean.fa"
    OUT="$OUT_DIR/${SAMPLE}.coverm.clean.tsv"

    log "----------------------------------------"
    log "Processing: $SAMPLE"

    [[ -s "$COV" ]] || {
        log "Missing CoverM file → skip"
        return
    }

    [[ -s "$CLEAN_CONTIGS" ]] || {
        log "Missing cleaned contigs → skip"
        return
    }

    # -------------------- SKIP --------------------
    if [[ -s "$OUT" && "${FORCE:-0}" != "1" ]]; then
        log "Clean CoverM exists → skipping"
        return
    fi

    TMP="$OUT_DIR/${SAMPLE}.tmp.tsv"

    # -------------------- STEP 1: BASIC CLEAN --------------------
    log "Filtering valid coverage entries"

    awk '
    NR==1 {next} 
    $3 > 0 && $3 != "NA" {
        print $1"\t"$3
    }
    ' "$COV" > "$TMP"

    # -------------------- STEP 2: MATCH CONTIG NAMES --------------------
    log "Matching with cleaned contigs"

    grep "^>" "$CLEAN_CONTIGS" | sed 's/^>//' | sort > "$OUT_DIR/${SAMPLE}.contigs.list"
    cut -f1 "$TMP" | sort > "$OUT_DIR/${SAMPLE}.coverm.list"

    comm -12 \
        "$OUT_DIR/${SAMPLE}.contigs.list" \
        "$OUT_DIR/${SAMPLE}.coverm.list" \
        > "$OUT_DIR/${SAMPLE}.shared.list"

    N_SHARED=$(wc -l < "$OUT_DIR/${SAMPLE}.shared.list")

    if [[ "$N_SHARED" -lt 10 ]]; then
        log "ERROR: Too few matching contigs ($N_SHARED)"
        exit 1
    fi

    log "Matching contigs: $N_SHARED"

    # -------------------- STEP 3: FILTER FINAL --------------------
    grep -F -w -f "$OUT_DIR/${SAMPLE}.shared.list" "$TMP" > "$OUT"

    [[ -s "$OUT" ]] || {
        log "ERROR: Clean CoverM empty"
        exit 1
    }

    rm -f "$TMP"

    log "Final entries: $(wc -l < "$OUT")"
    log "Finished: $SAMPLE"
}

# EXECUTION

if [[ "$process_sample" == true ]]; then
    run_clean_coverm "$COV" "$SAMPLE"
else
    for COV in "${COV_FILES[@]}"; do
        SAMPLE=$(basename "$COV" .coverm.tsv)
        run_clean_coverm "$COV" "$SAMPLE"
    done
fi

log "CoverM cleaning completed"
;;

# ============================================================
### BINNING (MAXBIN2 — CLEAN ASSEMBLY + CLEAN COVERM)
# ============================================================
binning_maxbin_prok)

set -euo pipefail
shopt -s nullglob

log "Starting MaxBin2 (clean assembly + clean CoverM)"

THREADS="${THREADS:-8}"
MIN_LEN="${MIN_LEN:-1500}"

ASM_DIR="$OUTPUT_DIR/assembly_clean"
COV_DIR="$OUTPUT_DIR/coverage_coverm_clean"

BIN_DIR="$OUTPUT_DIR/binning_maxbin_prok"
LOG_DIR="$OUTPUT_DIR/logs/binning_maxbin_prok"

mkdir -p "$BIN_DIR" "$LOG_DIR"

# -------------------- COLLECT SAMPLES --------------------
mapfile -t ASM_FILES < <(
    for d in "$ASM_DIR"/*; do
        [[ -d "$d" ]] || continue
        s=$(basename "$d")
        f="$d/${s}.contigs.clean.fa"
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total assemblies: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No cleaned assemblies found"; exit 1; }

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID > N )); then
    log "No sample → exiting"
    exit 0
fi

# -------------------- SELECT SAMPLE --------------------
CONTIGS="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$CONTIGS" .contigs.clean.fa)

COVERM="$COV_DIR/${SAMPLE}.coverm.clean.tsv"
ABUND="$BIN_DIR/${SAMPLE}/${SAMPLE}.maxbin.tsv"

OUTDIR="$BIN_DIR/$SAMPLE"

log "----------------------------------------"
log "Processing: $SAMPLE"

# -------------------- INPUT CHECK --------------------
[[ -s "$CONTIGS" ]] || { log "Missing contigs → skip"; exit 0; }
[[ -s "$COVERM" ]]  || { log "Missing clean coverm → skip"; exit 0; }

# HANDLE INCOMPLETE OUTPUT

if [[ -d "$OUTDIR" ]] && ! compgen -G "$OUTDIR/bin.*.fasta" > /dev/null; then
    log "Incomplete run detected → cleaning"
    rm -rf "$OUTDIR"
fi

# -------------------- SKIP --------------------
if [[ -d "$OUTDIR" ]] && compgen -G "$OUTDIR/bin.*.fasta" > /dev/null \
   && [[ "${FORCE:-0}" != "1" ]]; then
    log "MaxBin2 already done → skipping"
    exit 0
fi

# -------------------- CLEAN --------------------
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"


# FILTER CONTIGS BY LENGTH

FILTERED_CONTIGS="$OUTDIR/${SAMPLE}.min${MIN_LEN}.fa"

log "Filtering contigs (min length = $MIN_LEN)"

seqkit seq -m "$MIN_LEN" "$CONTIGS" > "$FILTERED_CONTIGS"

[[ -s "$FILTERED_CONTIGS" ]] || {
    log "ERROR: No contigs after filtering"
    exit 1
}

# PREPARE MAXBIN ABUNDANCE FILE (ALREADY CLEAN)

log "Preparing abundance file for MaxBin"

cp "$COVERM" "$ABUND"

[[ -s "$ABUND" ]] || {
    log "ERROR: Abundance file missing"
    exit 1
}

log "Abundance entries: $(wc -l < "$ABUND")"

# QUICK CONSISTENCY CHECK (SAFE, NON-BLOCKING)

log "Checking contig overlap (sanity check)"

grep "^>" "$FILTERED_CONTIGS" | sed 's/^>//' | sort > "$OUTDIR/contigs.list"
cut -f1 "$ABUND" | sort > "$OUTDIR/abund.list"

comm -12 "$OUTDIR/contigs.list" "$OUTDIR/abund.list" > "$OUTDIR/shared.list"

N_SHARED=$(wc -l < "$OUTDIR/shared.list")

if [[ "$N_SHARED" -lt 100 ]]; then
    log "WARNING: Low overlap ($N_SHARED) — check upstream steps"
else
    log "Matching contigs: $N_SHARED"
fi

# RUN MAXBIN2

log "Running MaxBin2"

run_MaxBin.pl \
    -contig "$FILTERED_CONTIGS" \
    -abund "$ABUND" \
    -out "$OUTDIR/bin" \
    -thread "$THREADS" \
    &>> "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
N_BINS=$(ls "$OUTDIR"/bin.*.fasta 2>/dev/null | wc -l || true)

if [[ "$N_BINS" -eq 0 ]]; then
    log "WARNING: No bins produced"
else
    log "Bins produced: $N_BINS"
fi

log "Finished MaxBin2 for $SAMPLE"
;;

# ============================================================
### BIN REFINEMENT (DAS TOOL — PROK, FINAL ROBUST VERSION)
# ============================================================
das_tool_prok)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting DAS Tool (prok bin refinement)"

THREADS="${THREADS:-16}"
MIN_LEN="${MIN_LEN:-1500}"

MB_DIR="$OUTPUT_DIR/binning_metabat_prok"
MX_DIR="$OUTPUT_DIR/binning_maxbin_prok"
ASM_CLEAN_DIR="$OUTPUT_DIR/assembly_clean"

DAS_DIR="$OUTPUT_DIR/das_tool_prok"
LOG_DIR="$OUTPUT_DIR/logs/das_tool_prok"

mkdir -p "$DAS_DIR" "$LOG_DIR"

# ============================================================
# DETECT INCOMPLETE SAMPLES (STRICT CHECK)
# ============================================================
mapfile -t SAMPLES < <(
    for d in "$MB_DIR"/*; do
        [[ -d "$d" ]] || continue

        SAMPLE=$(basename "$d")

        OUTDIR="$DAS_DIR/$SAMPLE"
        PREFIX="$OUTDIR/${SAMPLE}"

        SUMMARY="${PREFIX}_DASTool_summary.tsv"
        CONTIG2BIN="${PREFIX}_DASTool_contig2bin.tsv"
        BINS_DIR="${PREFIX}_DASTool_bins"

        COMPLETE=0

        if [[ -s "$SUMMARY" && -s "$CONTIG2BIN" && -d "$BINS_DIR" ]]; then

            N_BINS=$(find "$BINS_DIR" -name "*.fa" | wc -l)

            if [[ "$N_BINS" -gt 0 ]]; then
                COMPLETE=1
            fi
        fi

        if [[ "$COMPLETE" -eq 0 ]]; then
            echo "$SAMPLE"
        fi

    done | sort
)

N=${#SAMPLES[@]}

log "Samples needing DAS Tool: $N"

[[ "$N" -gt 0 ]] || {
    log "Nothing to run — all samples complete"
    exit 0
}

# ============================================================
# SINGLE-SAMPLE OVERRIDE
# ============================================================
if [[ -n "${TARGET_SAMPLE:-}" ]]; then

    log "TARGET_SAMPLE detected: $TARGET_SAMPLE"

    SAMPLES=("$TARGET_SAMPLE")
    N=1

    SLURM_ARRAY_TASK_ID=1
fi

# ============================================================
# ARRAY CHECK
# ============================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "No sample for this index → exiting"
    exit 0
fi

# ============================================================
# SELECT SAMPLE
# ============================================================
SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

MB="$MB_DIR/$SAMPLE"
MX="$MX_DIR/$SAMPLE"

ASM_CLEAN="$ASM_CLEAN_DIR/$SAMPLE/${SAMPLE}.contigs.clean.fa"

OUTDIR="$DAS_DIR/$SAMPLE"
PREFIX="$OUTDIR/${SAMPLE}"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# ============================================================
# INPUT CHECK
# ============================================================
[[ -d "$MB" && -d "$MX" ]] || {
    log "Missing MetaBAT or MaxBin directory → skipping"
    exit 0
}

[[ -s "$ASM_CLEAN" ]] || {
    log "Missing clean contigs file → skipping"
    exit 0
}

# ============================================================
# FINAL OUTPUT CHECK (EXTRA SAFETY)
# ============================================================
SUMMARY_FILE="${PREFIX}_DASTool_summary.tsv"
CONTIG2BIN_FILE="${PREFIX}_DASTool_contig2bin.tsv"
BINS_DIR="${PREFIX}_DASTool_bins"

if [[ -s "$SUMMARY_FILE" && -s "$CONTIG2BIN_FILE" && -d "$BINS_DIR" ]]; then

    N_BINS=$(find "$BINS_DIR" -name "*.fa" | wc -l)

    if [[ "$N_BINS" -gt 0 ]]; then
        log "Sample already completed successfully → skipping"
        exit 0
    fi
fi

# ============================================================
# CLEAN PARTIAL OUTPUTS ONLY
# ============================================================
if [[ -d "$OUTDIR" ]]; then
    log "Removing incomplete/partial output directory"
    rm -rf "$OUTDIR"
fi

mkdir -p "$OUTDIR"

# ============================================================
# FILTER CONTIGS (MATCH BINNING LENGTH)
# ============================================================
FILTERED_CONTIGS="$OUTDIR/${SAMPLE}.min${MIN_LEN}.fa"

log "Filtering contigs (minimum length = $MIN_LEN bp)"

seqkit seq -m "$MIN_LEN" "$ASM_CLEAN" > "$FILTERED_CONTIGS"

[[ -s "$FILTERED_CONTIGS" ]] || {
    log "ERROR: No contigs remaining after filtering"
    exit 1
}

# ============================================================
# CREATE SCAFFOLD2BIN FILES
# ============================================================
log "Creating scaffold2bin mapping files"

MB_MAP="$OUTDIR/metabat_scaffolds2bin.tsv"
MX_MAP="$OUTDIR/maxbin_scaffolds2bin.tsv"

# -------------------- MetaBAT --------------------
for f in "$MB"/bin*.fa; do
    [[ -s "$f" ]] || continue

    BIN=$(basename "$f" .fa)

    grep "^>" "$f" | \
        sed 's/^>//' | \
        awk -v bin="$BIN" '{print $1"\t"bin}'

done | sort -u > "$MB_MAP"

# -------------------- MaxBin --------------------
for f in "$MX"/bin*.fasta; do
    [[ -s "$f" ]] || continue

    BIN=$(basename "$f" .fasta)

    grep "^>" "$f" | \
        sed 's/^>//' | \
        awk -v bin="$BIN" '{print $1"\t"bin}'

done | sort -u > "$MX_MAP"

# ============================================================
# CHECK INPUT MAPS
# ============================================================
if [[ ! -s "$MB_MAP" && ! -s "$MX_MAP" ]]; then
    log "ERROR: No scaffold2bin mappings generated"
    exit 1
fi

if [[ ! -s "$MB_MAP" || ! -s "$MX_MAP" ]]; then
    log "WARNING: One binning method missing"
fi

# ============================================================
# CONSISTENCY CHECK
# ============================================================
log "Checking contig consistency"

grep "^>" "$FILTERED_CONTIGS" | \
    sed 's/^>//' | \
    sort > "$OUTDIR/contigs.list"

cat "$MB_MAP" "$MX_MAP" | \
    cut -f1 | \
    sort -u > "$OUTDIR/bin_contigs.list"

comm -12 \
    "$OUTDIR/contigs.list" \
    "$OUTDIR/bin_contigs.list" \
    > "$OUTDIR/shared.list"

N_SHARED=$(wc -l < "$OUTDIR/shared.list")

log "Shared contigs between bins and assembly: $N_SHARED"

if [[ "$N_SHARED" -lt 100 ]]; then
    log "ERROR: Low overlap between binned contigs and filtered assembly"
    exit 1
fi

# ============================================================
# RUN DAS TOOL
# ============================================================
log "Running DAS Tool"

set +e

DAS_Tool \
    -i "$MB_MAP,$MX_MAP" \
    -c "$FILTERED_CONTIGS" \
    -o "$PREFIX" \
    --labels metabat,maxbin \
    --threads "$THREADS" \
    --search_engine diamond \
    --write_bins \
    --score_threshold 0.5 \
    &>> "$LOG_DIR/${SAMPLE}.log"

EXIT_CODE=$?

set -e

log "DAS Tool exit code: $EXIT_CODE"

if [[ "$EXIT_CODE" -ne 0 ]]; then
    log "DAS Tool crashed — showing last 50 log lines"

    tail -n 50 "$LOG_DIR/${SAMPLE}.log"

    exit 1
fi

# ============================================================
# VALIDATION
# ============================================================
SUMMARY_FILE="${PREFIX}_DASTool_summary.tsv"
CONTIG2BIN_FILE="${PREFIX}_DASTool_contig2bin.tsv"
BINS_DIR="${PREFIX}_DASTool_bins"

if [[ -s "$SUMMARY_FILE" && -s "$CONTIG2BIN_FILE" && -d "$BINS_DIR" ]]; then

    N_BINS=$(find "$BINS_DIR" -name "*.fa" | wc -l)

    if [[ "$N_BINS" -gt 0 ]]; then
        log "SUCCESS: DAS Tool produced $N_BINS refined bins"
    else
        log "ERROR: No refined bins produced"
        exit 1
    fi

else
    log "ERROR: DAS Tool output incomplete"
    exit 1
fi

log "Finished DAS Tool for $SAMPLE"

;;

# ============================================================
### CHECKM (PROK QUALITY — FROM DAS TOOL, PER-SAMPLE)
# ============================================================
checkm_prok)

set -euo pipefail
shopt -s nullglob

log "Starting CheckM (from DAS Tool)"

THREADS="${THREADS:-16}"

DAS_DIR="$OUTPUT_DIR/das_tool_prok"
CHECKM_DIR="$OUTPUT_DIR/checkm_prok"
LOG_DIR="$OUTPUT_DIR/logs/checkm_prok"

mkdir -p "$CHECKM_DIR" "$LOG_DIR"

# ============================================================
# COLLECT VALID DAS TOOL SAMPLES
# ============================================================
mapfile -t SAMPLES < <(

    for d in "$DAS_DIR"/*; do

        [[ -d "$d" ]] || continue

        SAMPLE=$(basename "$d")

        BINS_DIR="$d/${SAMPLE}_DASTool_bins"

        if [[ -d "$BINS_DIR" ]] && \
           compgen -G "$BINS_DIR/*.fa" > /dev/null; then

            echo "$SAMPLE"
        fi

    done | sort
)

# ============================================================
# COUNT SAMPLES
# ============================================================
N=${#SAMPLES[@]}

log "Total samples detected: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No valid DAS Tool samples found"
    exit 1
}

# ============================================================
# ARRAY CHECK
# ============================================================
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
    log "ERROR: Run as SLURM array job"
    exit 1
}

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "No sample assigned to this array index → exiting"
    exit 0
fi

# ============================================================
# SELECT SAMPLE
# ============================================================
SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

BIN_DIR="$DAS_DIR/$SAMPLE/${SAMPLE}_DASTool_bins"

OUTDIR="$CHECKM_DIR/$SAMPLE"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# ============================================================
# INPUT CHECK
# ============================================================
[[ -d "$BIN_DIR" ]] || {
    log "Missing DAS Tool bins directory → skipping"
    exit 0
}

BINS=("$BIN_DIR"/*.fa)

if ! compgen -G "$BIN_DIR/*.fa" > /dev/null; then
    log "No DAS Tool bins found → skipping"
    exit 0
fi

# ============================================================
# SKIP COMPLETED
# ============================================================
if [[ -s "$OUTDIR/storage/bin_stats_ext.tsv" &&
      -s "$OUTDIR/${SAMPLE}_qa.tsv" &&
      "${FORCE:-0}" != "1" ]]; then

    log "CheckM already completed successfully → skipping"
    exit 0
fi

# ============================================================
# CLEAN PARTIAL OUTPUTS
# ============================================================
if [[ -d "$OUTDIR" ]]; then
    log "Removing incomplete CheckM output"
    rm -rf "$OUTDIR"
fi

mkdir -p "$OUTDIR/input"

# ============================================================
# FILTER SMALL BINS
# ============================================================
MIN_SIZE=50000

N_KEPT=0

log "Filtering bins smaller than ${MIN_SIZE} bp"

for BIN in "${BINS[@]}"; do

    SIZE=$(grep -v "^>" "$BIN" | tr -d '\n' | wc -c)

    if (( SIZE >= MIN_SIZE )); then
        cp "$BIN" "$OUTDIR/input/"
        ((++N_KEPT))
    fi
done

log "Bins retained after filtering: $N_KEPT"

[[ "$N_KEPT" -gt 0 ]] || {
    log "No bins passed minimum size filter → skipping"
    exit 0
}

# ============================================================
# RUN CHECKM
# ============================================================
log "Running CheckM lineage_wf"
log "This may take considerable time"

set +e

checkm lineage_wf \
    -x fa \
    "$OUTDIR/input" \
    "$OUTDIR" \
    -t "$THREADS" \
    &>> "$LOG_DIR/${SAMPLE}.log"

EXIT_CODE=$?

set -e

log "CheckM exit code: $EXIT_CODE"

if [[ "$EXIT_CODE" -ne 0 ]]; then

    log "CheckM failed — showing last 50 log lines"

    tail -n 50 "$LOG_DIR/${SAMPLE}.log"

    exit 1
fi

# ============================================================
# VALIDATE OUTPUTS
# ============================================================
[[ -s "$OUTDIR/storage/bin_stats_ext.tsv" &&
   -s "$OUTDIR/lineage.ms" ]] || {

    log "ERROR: CheckM outputs incomplete"

    exit 1
}

log "Finished CheckM for $SAMPLE"

;;

# ============================================================
### CHECKM QA REGENERATION (FAST RE-RUN)
# ============================================================
checkm_qa_prok)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting CheckM QA generation"

CHECKM_DIR="$OUTPUT_DIR/checkm_prok"

# ============================================================
# LOAD MODULES
# ============================================================
module --force purge
module load cluster/genius/login
module load CheckM/1.2.3-foss-2023b

# ============================================================
# TOOL CHECK
# ============================================================
command -v checkm >/dev/null || {
    log "ERROR: checkm not found"
    exit 1
}

# ============================================================
# INPUT CHECK
# ============================================================
[[ -d "$CHECKM_DIR" ]] || {
    log "ERROR: Missing CheckM directory"
    exit 1
}

# ============================================================
# LOOP OVER SAMPLES
# ============================================================
for OUTDIR in "$CHECKM_DIR"/*; do

    [[ -d "$OUTDIR" ]] || continue

    SAMPLE=$(basename "$OUTDIR")

    LINEAGE="$OUTDIR/lineage.ms"

    [[ -s "$LINEAGE" ]] || {
        log "Missing lineage.ms for $SAMPLE → skipping"
        continue
    }

    [[ -s "$OUTDIR/storage/bin_stats_ext.tsv" ]] || {
        log "Missing bin_stats_ext.tsv for $SAMPLE → skipping"
        continue
    }

    log "----------------------------------------"
    log "Generating QA table for: $SAMPLE"

    checkm qa \
        "$LINEAGE" \
        "$OUTDIR" \
        -o 1 \
        > "$OUTDIR/${SAMPLE}_qa.tsv"

    [[ -s "$OUTDIR/${SAMPLE}_qa.tsv" ]] || {
        log "ERROR: QA generation failed for $SAMPLE"
        exit 1
    }

    log "Finished QA generation for $SAMPLE"

done

log "Finished CheckM QA regeneration"

;;

# ============================================================
### CHECKM CLASSIFICATION (FROM DAS TOOL BINS)
# ============================================================
checkm_classify_prok)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

log "Starting CheckM classification"

CHECKM_DIR="$OUTPUT_DIR/checkm_prok"
DAS_DIR="$OUTPUT_DIR/das_tool_prok"

CLASS_DIR="$OUTPUT_DIR/checkm_classified_prok"

SUMMARY="$CLASS_DIR/checkm_summary.tsv"

# ============================================================
# PREP OUTPUT DIRECTORIES
# ============================================================
mkdir -p "$CLASS_DIR"/{HQ,MQ,LQ}

# Clean previous outputs
rm -f "$SUMMARY"

rm -rf "$CLASS_DIR"/HQ/*
rm -rf "$CLASS_DIR"/MQ/*
rm -rf "$CLASS_DIR"/LQ/*

# ============================================================
# HEADER
# ============================================================
echo -e "Sample\tBin\tCompleteness\tContamination\tCategory" \
    > "$SUMMARY"

# ============================================================
# LOOP OVER SAMPLES
# ============================================================
for SAMPLE_DIR in "$CHECKM_DIR"/*; do

    [[ -d "$SAMPLE_DIR" ]] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")

    QA_FILE="$SAMPLE_DIR/${SAMPLE}_qa.tsv"

    BIN_DIR="$DAS_DIR/$SAMPLE/${SAMPLE}_DASTool_bins"

    # --------------------------------------------------------
    # INPUT CHECKS
    # --------------------------------------------------------
    [[ -s "$QA_FILE" ]] || {
        log "Missing QA file for $SAMPLE → skipping"
        continue
    }

    [[ -d "$BIN_DIR" ]] || {
        log "Missing DAS Tool bins for $SAMPLE → skipping"
        continue
    }

    log "Processing sample: $SAMPLE"

    TMPFILE="$CLASS_DIR/tmp_${SAMPLE}.tsv"

    rm -f "$TMPFILE"

    # ========================================================
    # PARSE CHECKM QA FILE
    # ========================================================
    awk \
    -v sample="$SAMPLE" \
    -v out="$SUMMARY" \
    -v tmp="$TMPFILE" \
    '
    BEGIN {
        OFS="\t"
    }

    # --------------------------------------------------------
    # Skip logs / separators / headers
    # --------------------------------------------------------
    /^\[/ { next }
    /^---/ { next }
    /Bin Id/ { next }

    {
        line = $0

        # Remove leading whitespace
        sub(/^[[:space:]]+/, "", line)

        # Split line
        n = split(line, a, /[[:space:]]+/)

        if (n < 5)
            next

        bin = a[1]

        # ----------------------------------------------------
        # Keep only real bins
        # ----------------------------------------------------
        if (bin !~ /bin/)
            next

        # ----------------------------------------------------
        # Extract numeric columns
        # ----------------------------------------------------
        nnums = 0

        for (i = 1; i <= n; i++) {

            if (a[i] ~ /^[0-9.]+$/) {

                nums[++nnums] = a[i]
            }
        }

        # Need:
        # completeness contamination strain_heterogeneity
        if (nnums < 3)
            next

        # ----------------------------------------------------
        # Final numeric columns:
        # completeness contamination strain_heterogeneity
        # ----------------------------------------------------
        comp = nums[nnums-2] + 0
        cont = nums[nnums-1] + 0

        # ----------------------------------------------------
        # Classification
        # ----------------------------------------------------
        if (comp >= 90 && cont <= 5)
            cat = "HQ"
        else if (comp >= 50 && cont <= 10)
            cat = "MQ"
        else
            cat = "LQ"

        # ----------------------------------------------------
        # Write summary
        # ----------------------------------------------------
        print sample, bin, comp, cont, cat >> out

        # Temporary mapping for copying bins
        print bin, cat >> tmp
    }
    ' "$QA_FILE"

    # ========================================================
    # COPY CLASSIFIED BINS
    # ========================================================
    if [[ -s "$TMPFILE" ]]; then

        while IFS=$'\t' read -r BIN CATEGORY; do

            MATCH=$(find "$BIN_DIR" -name "${BIN}.fa" | head -n 1)

            if [[ -n "$MATCH" && -s "$MATCH" ]]; then

                cp "$MATCH" \
                   "$CLASS_DIR/$CATEGORY/${SAMPLE}_$(basename "$MATCH")"

            else
                log "WARNING: Missing source bin for $BIN"
            fi

        done < "$TMPFILE"

        rm -f "$TMPFILE"
    fi

done

# ============================================================
# FINAL COUNTS
# ============================================================
HQ_N=$(find "$CLASS_DIR/HQ" -name "*.fa" | wc -l)
MQ_N=$(find "$CLASS_DIR/MQ" -name "*.fa" | wc -l)
LQ_N=$(find "$CLASS_DIR/LQ" -name "*.fa" | wc -l)

log "HQ bins: $HQ_N"
log "MQ bins: $MQ_N"
log "LQ bins: $LQ_N"

log "Finished CheckM classification"

;;

# ============================================================
### CHECKM SUMMARY PER SAMPLE
# ============================================================
checkm_summary_sample_prok)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

IN="$OUTPUT_DIR/checkm_classified_prok/checkm_summary.tsv"

CHECKM_DIR="$OUTPUT_DIR/checkm_prok"

OUT="$OUTPUT_DIR/checkm_classified_prok/checkm_summary_per_sample.tsv"

log "Creating CheckM per-sample summary"

# ============================================================
# GET ALL SAMPLES
# ============================================================
mapfile -t ALL_SAMPLES < <(

    for d in "$CHECKM_DIR"/*; do

        [[ -d "$d" ]] || continue

        basename "$d"

    done | sort
)

# ============================================================
# BUILD SUMMARY
# ============================================================
awk '
BEGIN {
    FS="\t"
    OFS="\t"
}

NR == 1 { next }

{
    sample = $1
    comp_v = $3 + 0
    cont_v = $4 + 0
    cat = $5

    comp[sample] += comp_v
    cont[sample] += cont_v
    n[sample]++

    if (cat == "HQ")
        HQ[sample]++
    else if (cat == "MQ")
        MQ[sample]++
    else if (cat == "LQ")
        LQ[sample]++
}

END {

    for (s in n) {

        mean_comp[s] = comp[s] / n[s]
        mean_cont[s] = cont[s] / n[s]
    }

    for (s in mean_comp) {

        print s,
              mean_comp[s],
              mean_cont[s],
              n[s],
              (s in HQ ? HQ[s] : 0),
              (s in MQ ? MQ[s] : 0),
              (s in LQ ? LQ[s] : 0)
    }
}
' "$IN" > "$OUT.tmp"

# ============================================================
# FINAL TABLE WITH ALL SAMPLES
# ============================================================
{
    echo -e "sample\tmean_completeness\tmean_contamination\tn_bins\tHQ\tMQ\tLQ"

    for SAMPLE in "${ALL_SAMPLES[@]}"; do

        LINE=$(awk -F'\t' -v s="$SAMPLE" '$1==s' "$OUT.tmp")

        if [[ -n "$LINE" ]]; then

            echo "$LINE"

        else

            echo -e "${SAMPLE}\t0\t0\t0\t0\t0\t0"
        fi
    done

} > "$OUT"

rm -f "$OUT.tmp"

log "Created: $OUT"

;;

# ============================================================
### BAT (BIN-LEVEL — PROKARYOTIC MAGS)
# ============================================================
bat_prok)

set -euo pipefail
shopt -s nullglob

log() {
    echo "[$(date '+%F %T')] $*"
}

log "Starting BAT (prok MAG taxonomy)"

# ============================================================
# SETTINGS
# ============================================================
THREADS="${THREADS:-8}"

CLASS_DIR="$OUTPUT_DIR/checkm_classified_prok"

BAT_DIR="$OUTPUT_DIR/bat_prok"

LOG_DIR="$OUTPUT_DIR/logs/bat_prok"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$BAT_DIR" "$LOG_DIR"

# ============================================================
# DEPENDENCY CHECKS
# ============================================================
command -v "$CAT_CMD" >/dev/null || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null || {
    log "ERROR: DIAMOND not found"
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
    log "ERROR: Missing CAT database"
    exit 1
}

[[ -d "$CAT_TAX" ]] || {
    log "ERROR: Missing CAT taxonomy"
    exit 1
}

# ============================================================
# COLLECT SAMPLES
# ============================================================
mapfile -t SAMPLES < <(

    find "$CLASS_DIR"/HQ \
         "$CLASS_DIR"/MQ \
         "$CLASS_DIR"/LQ \
         -maxdepth 1 \
         -type f \
         -name "*.fa" \
         2>/dev/null \
    | xargs -n 1 basename \
    | cut -d'_' -f1 \
    | sort -u
)

N=${#SAMPLES[@]}

log "Total samples detected: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No classified MAG samples found"
    exit 1
}

# ============================================================
# ARRAY CHECK
# ============================================================
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
    log "ERROR: Run as SLURM array job"
    exit 1
}

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then

    log "No sample assigned to this array index → exiting"

    exit 0
fi

# ============================================================
# SELECT SAMPLE
# ============================================================
SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

OUTDIR="$BAT_DIR/$SAMPLE"

PREFIX="$OUTDIR/${SAMPLE}"

BINS_TMP="$OUTDIR/bins"

# ============================================================
# LOCAL TMPDIR (VERY IMPORTANT FOR DIAMOND)
# ============================================================
SAMPLE_TMP="${TMPDIR:-/tmp}/${SAMPLE}"

mkdir -p "$SAMPLE_TMP"

export TMPDIR="$SAMPLE_TMP"
export TEMP="$SAMPLE_TMP"
export TMP="$SAMPLE_TMP"

log "Using sample TMPDIR: $TMPDIR"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# ============================================================
# SKIP COMPLETED
# ============================================================
if [[ -s "${PREFIX}.bin2classification.txt" &&
      -s "${PREFIX}.ORF2LCA.txt" &&
      "${FORCE:-0}" != "1" ]]; then

    log "BAT already completed successfully → skipping"

    exit 0
fi

# ============================================================
# CLEAN PARTIAL OUTPUTS
# ============================================================
if [[ -d "$OUTDIR" ]]; then

    log "Removing incomplete BAT output"

    rm -rf "$OUTDIR"
fi

mkdir -p "$BINS_TMP"

# ============================================================
# COLLECT SAMPLE MAGS
# ============================================================
N_BINS=0

mapfile -t SAMPLE_BINS < <(

    find "$CLASS_DIR"/HQ \
         "$CLASS_DIR"/MQ \
         "$CLASS_DIR"/LQ \
         -maxdepth 1 \
         -type f \
         -name "${SAMPLE}_*.fa" \
         | sort
)

for BIN in "${SAMPLE_BINS[@]}"; do

    [[ -s "$BIN" ]] || continue

    cp "$BIN" "$BINS_TMP/"

    ((++N_BINS))
done

log "Bins copied for sample: $N_BINS"

[[ "$N_BINS" -gt 0 ]] || {
    log "ERROR: No MAGs found for sample"
    exit 1
}

# ============================================================
# DEBUG
# ============================================================
log "Files inside BAT bins directory:"

ls -lh "$BINS_TMP"

# ============================================================
# RUN BAT
# ============================================================
log "Running CAT_pack BAT"

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
    2>&1 | tee "$LOG_DIR/${SAMPLE}.log"

EXIT_CODE=${PIPESTATUS[0]}

set -e

log "BAT exit code: $EXIT_CODE"

# ============================================================
# HANDLE FAILURE
# ============================================================
if [[ "$EXIT_CODE" -ne 0 ]]; then

    log "BAT failed — showing last 50 log lines"

    tail -n 50 "$LOG_DIR/${SAMPLE}.log"

    exit 1
fi

# ============================================================
# VALIDATE OUTPUTS
# ============================================================
CLASS_FILE="${PREFIX}.bin2classification.txt"

ORF_FILE="${PREFIX}.ORF2LCA.txt"

[[ -s "$CLASS_FILE" ]] || {

    log "ERROR: Missing bin2classification output"

    exit 1
}

[[ -s "$ORF_FILE" ]] || {

    log "ERROR: Missing ORF2LCA output"

    exit 1
}

# ============================================================
# CLEAN TMP FILES
# ============================================================
rm -rf "$TMPDIR"

# ============================================================
# FINAL MESSAGE
# ============================================================
log "Finished BAT for sample: $SAMPLE"

;;

# ============================================================
### GTDB-TK (PROK TAXONOMY — PER SAMPLE, ARRAY-AWARE)
# ============================================================
gtdbtk)

set -euo pipefail
shopt -s nullglob

log "Starting GTDB-Tk classification"

THREADS=8
PPLACER_CPUS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CLASS_DIR="$OUTPUT_DIR/checkm_classified_prok"

GTDB_DIR="$OUTPUT_DIR/gtdbtk"

LOG_DIR="$OUTPUT_DIR/logs/gtdbtk"

# ============================================================
# GTDB DATABASE
# ============================================================
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-/staging/leuven/stg_00100/Rodgee/databases/gtdbtk/release232}"

log "GTDBTK_DATA_PATH: $GTDBTK_DATA_PATH"

[[ -d "$GTDBTK_DATA_PATH" ]] || {
    log "ERROR: GTDBTK_DATA_PATH not found"
    exit 1
}

mkdir -p "$GTDB_DIR" "$LOG_DIR"

# ============================================================
# CHECK GTDB-TK
# ============================================================
command -v gtdbtk >/dev/null || {
    log "ERROR: gtdbtk not found"
    exit 1
}

# ============================================================
# COLLECT SAMPLES
# ============================================================
mapfile -t SAMPLES < <(

    find "$CLASS_DIR"/HQ \
         "$CLASS_DIR"/MQ \
         "$CLASS_DIR"/LQ \
         -maxdepth 1 \
         -type f \
         -name "*.fa" \
         2>/dev/null \
    | xargs -n 1 basename \
    | cut -d'_' -f1 \
    | sort -u
)

N=${#SAMPLES[@]}

log "Total samples detected: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No samples detected"
    exit 1
}

# ============================================================
# ARRAY CHECK
# ============================================================
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
    log "ERROR: Run as SLURM array job"
    exit 1
}

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then

    log "No sample assigned to this array index → exiting"

    exit 0
fi

# ============================================================
# SELECT SAMPLE
# ============================================================
SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

OUTDIR="$GTDB_DIR/$SAMPLE"

TMP_INPUT="$OUTDIR/input"

log "----------------------------------------"
log "Processing sample: $SAMPLE"

# ============================================================
# OUTPUT FILES
# ============================================================
BAC="$OUTDIR/classify/gtdbtk.bac120.summary.tsv"

ARC="$OUTDIR/classify/gtdbtk.ar53.summary.tsv"

# ============================================================
# SKIP COMPLETED
# ============================================================
if [[ (
        -s "$BAC" || \
        -s "$ARC"
      ) && "${FORCE:-0}" != "1" ]]; then

    log "GTDB-Tk already completed successfully → skipping"

    exit 0
fi

# ============================================================
# CLEAN PARTIAL OUTPUTS
# ============================================================
if [[ -d "$OUTDIR" ]]; then

    log "Removing incomplete GTDB-Tk output"

    rm -rf "$OUTDIR"
fi

mkdir -p "$TMP_INPUT"

# ============================================================
# COLLECT SAMPLE MAGS
# ============================================================
N_MAGS=0

mapfile -t SAMPLE_MAGS < <(

    find "$CLASS_DIR"/HQ \
         "$CLASS_DIR"/MQ \
         "$CLASS_DIR"/LQ \
         -maxdepth 1 \
         -type f \
         -name "${SAMPLE}_*.fa" \
         | sort
)

for MAG in "${SAMPLE_MAGS[@]}"; do

    [[ -s "$MAG" ]] || continue

    BASE=$(basename "$MAG" .fa)

    # --------------------------------------------------------
    # DETERMINE QUALITY
    # --------------------------------------------------------
    if [[ "$MAG" == *"/HQ/"* ]]; then
        QUALITY="HQ"

    elif [[ "$MAG" == *"/MQ/"* ]]; then
        QUALITY="MQ"

    else
        QUALITY="LQ"
    fi

    NEW_NAME="${BASE}_${QUALITY}.fa"

    cp "$MAG" "$TMP_INPUT/$NEW_NAME"

    ((++N_MAGS))
done

log "MAGs copied for sample: $N_MAGS"

[[ "$N_MAGS" -gt 0 ]] || {
    log "ERROR: No MAGs found for sample"
    exit 1
}

# ============================================================
# DEBUG
# ============================================================
log "Files inside TMP_INPUT:"
ls -lh "$TMP_INPUT"

log "THREADS: $THREADS"
log "PPLACER_CPUS: $PPLACER_CPUS"

# ============================================================
# TMPDIR
# ============================================================
GTDB_TMP="${TMPDIR:-$OUTDIR/tmp}"

mkdir -p "$GTDB_TMP"

export TMPDIR="$GTDB_TMP"

log "TMPDIR: $TMPDIR"

# ============================================================
# RUN GTDB-TK
# ============================================================
log "Running GTDB-Tk classify_wf"

set +e

gtdbtk classify_wf \
    --genome_dir "$TMP_INPUT" \
    --out_dir "$OUTDIR" \
    --extension fa \
    --cpus "$THREADS" \
    --pplacer_cpus "$PPLACER_CPUS" \
    --tmpdir "$TMPDIR" \
    &>> "$LOG_DIR/${SAMPLE}.log"

EXIT_CODE=$?

set -e

log "GTDB-Tk exit code: $EXIT_CODE"

# ============================================================
# HANDLE FAILURE
# ============================================================
if [[ "$EXIT_CODE" -ne 0 ]]; then

    log "GTDB-Tk failed — showing last 50 log lines"

    tail -n 50 "$LOG_DIR/${SAMPLE}.log"

    exit 1
fi

# ============================================================
# VALIDATE OUTPUT
# ============================================================
if [[ ! -s "$BAC" && ! -s "$ARC" ]]; then

    log "ERROR: GTDB-Tk output incomplete"

    exit 1
fi

# ============================================================
# REPORT OUTPUTS
# ============================================================
if [[ -s "$BAC" ]]; then
    log "Bacterial classifications detected"
fi

if [[ -s "$ARC" ]]; then
    log "Archaeal classifications detected"
fi

# ============================================================
# FINAL MESSAGE
# ============================================================
log "Finished GTDB-Tk for sample: $SAMPLE"

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

    awk -v sample="$SAMPLE" '
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
    quality=$5

    # --------------------------------------------------------
    # REBUILD GTDB KEY
    # --------------------------------------------------------
    key=sample "_" bin "_" quality

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
### TIARA SPLIT (PER-SAMPLE, ARRAY-AWARE)
# ============================================================
tiara_split)

set -euo pipefail
shopt -s nullglob

log "Starting Tiara split step (array mode)"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-8}"
ASM_DIR="$OUTPUT_DIR/assembly"
TIARA_DIR="$OUTPUT_DIR/tiara"

mkdir -p "$TIARA_DIR"

# -------------------- COLLECT ASSEMBLIES --------------------
mapfile -t ASM_FILES < <(
    for d in "$ASM_DIR"/*; do
        [[ -d "$d" ]] || continue
        SAMPLE=$(basename "$d")
        ASM="$d/${SAMPLE}.contigs.fa"
        [[ -s "$ASM" ]] && echo "$ASM"
    done | sort
)

N=${#ASM_FILES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || {
    log "ERROR: No valid assemblies found"
    exit 1
}

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
CONTIGS="${ASM_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$CONTIGS" | sed 's/.contigs.fa//')

OUT="$TIARA_DIR/$SAMPLE"
mkdir -p "$OUT"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "Contigs: $CONTIGS"

# -------------------- CHECK TIARA --------------------
command -v tiara >/dev/null || {
    log "ERROR: tiara not found (activate conda in SLURM)"
    exit 1
}

tiara --help >/dev/null 2>&1 || {
    log "ERROR: tiara not working"
    exit 1
}

# -------------------- SKIP --------------------
if [[ -s "$OUT/tiara.txt" && "${FORCE:-0}" != "1" ]]; then
    log "Tiara already done → skipping"
    exit 0
fi

# -------------------- RUN TIARA --------------------
log "Running Tiara"

tiara \
    -i "$CONTIGS" \
    -o "$OUT/tiara.txt" \
    -t "$THREADS" \
    --min_len 3000

# -------------------- VALIDATE --------------------
[[ -s "$OUT/tiara.txt" ]] || {
    log "ERROR: Tiara failed"
    exit 1
}

# -------------------- SPLIT IDS --------------------
awk 'NR>1 && $5=="eukarya"{print $1}' "$OUT/tiara.txt" > "$OUT/euk_ids.txt"

awk 'NR>1 && ($5=="bacteria" || $5=="archaea" || $5=="prokarya"){print $1}' \
    "$OUT/tiara.txt" > "$OUT/prok_ids.txt"

awk 'NR>1 && $5=="organelle"{print $1}' "$OUT/tiara.txt" > "$OUT/organelle_ids.txt"

awk 'NR>1 && $5=="organelle" && $6=="mitochondrion"{print $1}' \
    "$OUT/tiara.txt" > "$OUT/mito_ids.txt"

awk 'NR>1 && $5=="organelle" && $6=="plastid"{print $1}' \
    "$OUT/tiara.txt" > "$OUT/plastid_ids.txt"

# -------------------- FASTA SPLIT (ONE-PASS, ROBUST) --------------------
log "Splitting FASTA (single-pass)"

awk -v euk="$OUT/euk_ids.txt" \
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
    while ((getline < euk)   > 0) euk_ids[$1]=1
    while ((getline < prok)  > 0) prok_ids[$1]=1
    while ((getline < org)   > 0) org_ids[$1]=1
    while ((getline < mito)  > 0) mito_ids[$1]=1
    while ((getline < plast) > 0) plast_ids[$1]=1

    close(euk); close(prok); close(org); close(mito); close(plast)
}

/^>/ {
    split(substr($0,2), a, " ")
    header = a[1]

    is_euk  = (header in euk_ids)
    is_prok = (header in prok_ids)
    is_org  = (header in org_ids)
    is_mito = (header in mito_ids)
    is_plast= (header in plast_ids)
}

{
    if (is_euk)  print >> euk_out
    if (is_prok) print >> prok_out
    if (is_org)  print >> org_out
    if (is_mito) print >> mito_out
    if (is_plast)print >> plast_out
}

END {
    close(euk_out)
    close(prok_out)
    close(org_out)
    close(mito_out)
    close(plast_out)
}
' "$CONTIGS"

# -------------------- EMPTY FILE SAFETY --------------------
for f in euk.fa prok.fa organelle.fa mito.fa plastid.fa; do
    [[ -s "$OUT/$f" ]] || : > "$OUT/$f"
done

# -------------------- SUMMARY --------------------
log "Summary for $SAMPLE:"
log "  euk:       $(wc -l < "$OUT/euk_ids.txt")"
log "  prok:      $(wc -l < "$OUT/prok_ids.txt")"
log "  organelle: $(wc -l < "$OUT/organelle_ids.txt")"
log "  mito:      $(wc -l < "$OUT/mito_ids.txt")"
log "  plastid:   $(wc -l < "$OUT/plastid_ids.txt")"

log "Finished Tiara split for $SAMPLE"
;;

# ============================================================
### CAT (CONTIG-LEVEL — TIARA PROKARYOTIC CONTIGS)
# ============================================================
cat_contigs_prok)

set -euo pipefail
shopt -s nullglob

log "Starting CAT (contig-level, prok)"

THREADS="${THREADS:-8}"

TARGET_SAMPLE="${TARGET_SAMPLE:-}"

TIARA_DIR="$OUTPUT_DIR/tiara"

CAT_DIR="$OUTPUT_DIR/cat_contigs_prok"

LOG_DIR="$OUTPUT_DIR/logs/cat_contigs_prok"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$CAT_DIR" "$LOG_DIR"

# ============================================================
# DEPENDENCY CHECK
# ============================================================
command -v "$CAT_CMD" >/dev/null || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null || {
    log "ERROR: DIAMOND not found"
    exit 1
}

# ============================================================
# FUNCTION TO PROCESS ONE SAMPLE
# ============================================================
run_cat() {

    local SAMPLE="$1"

    local SAMPLE_DIR="$TIARA_DIR/$SAMPLE"

    local CONTIGS="$SAMPLE_DIR/prok.fa"

    local OUTDIR="$CAT_DIR/$SAMPLE"

    local PREFIX="$OUTDIR/${SAMPLE}"

    log "----------------------------------------"
    log "Processing sample: $SAMPLE"

    # --------------------------------------------------------
    # CHECK INPUT
    # --------------------------------------------------------
    [[ -s "$CONTIGS" ]] || {

        log "No prok.fa found → skipping"

        return
    }

    # --------------------------------------------------------
    # SKIP COMPLETED
    # --------------------------------------------------------
    if [[ -s "${PREFIX}.contigs2classification.txt" &&
          "${FORCE:-0}" != "1" ]]; then

        log "CAT already completed → skipping"

        return
    fi

    # --------------------------------------------------------
    # CLEAN PARTIAL OUTPUTS
    # --------------------------------------------------------
    rm -rf "$OUTDIR"

    mkdir -p "$OUTDIR"

    # --------------------------------------------------------
    # DEBUG
    # --------------------------------------------------------
    log "Input prokaryotic contigs:"
    log "$CONTIGS"

    log "Output prefix:"
    log "$PREFIX"

    log "Threads: $THREADS"

    # --------------------------------------------------------
    # RUN CAT
    # --------------------------------------------------------
    log "Running CAT_pack contigs"

    set +e

    "$CAT_CMD" contigs \
        -c "$CONTIGS" \
        -d "$CAT_DB" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --nproc "$THREADS" \
        --block_size 2 \
        2>&1 | tee "$LOG_DIR/${SAMPLE}.log"

    EXIT_CODE=${PIPESTATUS[0]}

    set -e

    log "CAT exit code: $EXIT_CODE"

    # --------------------------------------------------------
    # HANDLE FAILURE
    # --------------------------------------------------------
    if [[ "$EXIT_CODE" -ne 0 ]]; then

        log "ERROR: CAT failed for $SAMPLE"

        tail -n 50 "$LOG_DIR/${SAMPLE}.log"

        return
    fi

    # --------------------------------------------------------
    # VALIDATE OUTPUT
    # --------------------------------------------------------
    [[ -s "${PREFIX}.contigs2classification.txt" ]] || {

        log "ERROR: Missing contigs2classification output"

        return
    }

    log "Finished CAT for sample: $SAMPLE"
}

# ============================================================
# MODE SELECTION
# ============================================================

# ------------------------------------------------------------
# MANUAL MODE
# ------------------------------------------------------------
if [[ -n "$TARGET_SAMPLE" ]]; then

    log "Manual mode → $TARGET_SAMPLE"

    run_cat "$TARGET_SAMPLE"

# ------------------------------------------------------------
# ARRAY MODE
# ------------------------------------------------------------
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    mapfile -t SAMPLES < <(

        for d in "$TIARA_DIR"/*; do

            [[ -d "$d" ]] || continue

            SAMPLE=$(basename "$d")

            CONTIGS="$d/prok.fa"

            [[ -s "$CONTIGS" ]] || continue

            PREFIX="$CAT_DIR/$SAMPLE/${SAMPLE}"

            # only incomplete samples
            [[ ! -s "${PREFIX}.contigs2classification.txt" ]] || continue

            echo "$SAMPLE"

        done | sort
    )

    N=${#SAMPLES[@]}

    log "Incomplete samples detected: $N"

    [[ "$N" -gt 0 ]] || {

        log "Nothing left to process"

        exit 0
    }

    if (( SLURM_ARRAY_TASK_ID > N )); then

        log "No sample assigned to this array index"

        exit 0
    fi

    SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"

    run_cat "$SAMPLE"

# ------------------------------------------------------------
# AUTO MODE
# ------------------------------------------------------------
else

    log "Auto mode → running only incomplete samples"

    for d in "$TIARA_DIR"/*; do

        [[ -d "$d" ]] || continue

        SAMPLE=$(basename "$d")

        CONTIGS="$d/prok.fa"

        PREFIX="$CAT_DIR/$SAMPLE/${SAMPLE}"

        if [[ -s "$CONTIGS" &&
              ! -s "${PREFIX}.contigs2classification.txt" ]]; then

            run_cat "$SAMPLE"

        else

            log "Skipping $SAMPLE (done or missing prok.fa)"
        fi
    done
fi

log "CAT prok contig step completed"

;;

# ============================================================
### REMAG (EUK BINNING — PER-SAMPLE, ARRAY-AWARE)
# ============================================================
remag)

set -euo pipefail
shopt -s nullglob

log "Starting reMAG step (array mode)"

# -------------------- CONFIG --------------------
THREADS="${THREADS:-8}"

TIARA_DIR="$OUTPUT_DIR/tiara"
MAP_DIR="$OUTPUT_DIR/mapping"
REMAG_DIR="$OUTPUT_DIR/remag"
LOG_DIR="$OUTPUT_DIR/logs/remag"

mkdir -p "$REMAG_DIR" "$LOG_DIR"

# -------------------- CHECK INPUT DIRS --------------------
[[ -d "$TIARA_DIR" ]] || { log "ERROR: Missing TIARA dir: $TIARA_DIR"; exit 1; }
[[ -d "$MAP_DIR"   ]] || { log "ERROR: Missing mapping dir: $MAP_DIR"; exit 1; }

# -------------------- TOOL CHECK --------------------
if ! command -v remag >/dev/null 2>&1; then
    log "ERROR: remag not found (activate conda in SLURM job)"
    exit 1
fi

# -------------------- COLLECT SAMPLES --------------------
mapfile -t TIARA_SAMPLES < <(
    for d in "$TIARA_DIR"/*; do
        [[ -d "$d" ]] && echo "$d"
    done | sort
)

N=${#TIARA_SAMPLES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No Tiara samples found"; exit 1; }

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT SAMPLE --------------------
COND="${TIARA_SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$COND")

EUK_FASTA="$COND/euk.fa"
BAM="$MAP_DIR/${SAMPLE}.bam"
OUTDIR="$REMAG_DIR/$SAMPLE"

log "----------------------------------------"
log "Processing sample: $SAMPLE"
log "EUK_FASTA: $EUK_FASTA"
log "BAM: $BAM"

# -------------------- INPUT VALIDATION --------------------
if [[ ! -s "$EUK_FASTA" ]]; then
    log "WARNING: No eukaryotic contigs → skipping"
    exit 0
fi

if [[ ! -s "$BAM" ]]; then
    log "WARNING: Missing BAM → skipping"
    exit 0
fi

# -------------------- SKIP --------------------
if [[ -s "$OUTDIR/bins.csv" && "${FORCE:-0}" != "1" ]]; then
    log "reMAG already done → skipping"
    exit 0
fi

# -------------------- CLEAN --------------------
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# -------------------- RUN REMAG --------------------
log "Running reMAG"

remag \
    "$EUK_FASTA" \
    -c "$BAM" \
    -o "$OUTDIR" \
    --skip-bacterial-filter \
    &>> "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
if [[ ! -s "$OUTDIR/bins.csv" ]]; then
    log "WARNING: No bins produced for $SAMPLE"
    exit 0
fi

# -------------------- SUMMARY --------------------
log "Finished reMAG for $SAMPLE"

if [[ -d "$OUTDIR/bins" ]]; then
    log "Number of bins: $(ls "$OUTDIR/bins"/*.fa 2>/dev/null | wc -l)"
else
    log "WARNING: No bins directory created"
fi

;;

# ============================================================
### BUSCO (EUK QUALITY — PER-BIN, ARRAY-AWARE)
# ============================================================
busco)

set -euo pipefail
shopt -s nullglob

log "Starting BUSCO step (array mode)"

THREADS="${THREADS:-8}"

REMAG_DIR="$OUTPUT_DIR/remag"
BUSCO_DIR="$OUTPUT_DIR/busco"
LOG_DIR="$OUTPUT_DIR/logs/busco"

mkdir -p "$BUSCO_DIR" "$LOG_DIR"

# -------------------- COLLECT ALL BINS --------------------
mapfile -t BIN_FILES < <(
    for f in "$REMAG_DIR"/*/bins/*.fa; do
        [[ -s "$f" ]] && echo "$f"
    done | sort
)

N=${#BIN_FILES[@]}

log "Total bins detected: $N"

if [[ "$N" -eq 0 ]]; then
    log "ERROR: No bins found"
    exit 1
fi

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
    log "ERROR: Task ID out of range (1–$N)"
    exit 1
fi

# -------------------- SELECT BIN --------------------
BIN="${BIN_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

SAMPLE=$(basename "$(dirname "$(dirname "$BIN")")")
BIN_NAME=$(basename "$BIN" .fa)

OUT_PREFIX="${SAMPLE}_${BIN_NAME}"
OUTDIR="$BUSCO_DIR/$OUT_PREFIX"

log "----------------------------------------"
log "Processing: $SAMPLE / $BIN_NAME"
log "BIN: $BIN"

# -------------------- FILTER SMALL BINS --------------------
SIZE=$(wc -c < "$BIN")
MIN_SIZE=50000  # 50 kb threshold (adjust if needed)

if (( SIZE < MIN_SIZE )); then
    log "WARNING: Bin too small (${SIZE} bp) → skipping BUSCO"
    exit 0
fi

# -------------------- SKIP --------------------
if compgen -G "$OUTDIR/short_summary.*.txt" > /dev/null && [[ "${FORCE:-0}" != "1" ]]; then
    log "BUSCO already done → skipping"
    exit 0
fi

# -------------------- CLEAN --------------------
rm -rf "$OUTDIR"

# -------------------- RUN BUSCO --------------------
log "Running BUSCO"

busco \
    -i "$BIN" \
    -l eukaryota_odb10 \
    -o "$OUT_PREFIX" \
    -m genome \
    -c "$THREADS" \
    --out_path "$BUSCO_DIR" \
    &>> "$LOG_DIR/${OUT_PREFIX}.log"

# -------------------- VALIDATE --------------------
if ! compgen -G "$OUTDIR/short_summary.*.txt" > /dev/null; then
    log "ERROR: BUSCO failed for $OUT_PREFIX"
    exit 1
fi

log "Finished BUSCO: $OUT_PREFIX"

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
### RAT (READ-LEVEL — TIARA EUK CONTIGS + EXISTING CAT OUTPUTS)
### RUN NEW SAMPLES + RERUN FAILED ONLY
# ============================================================
rat_reads)

set -euo pipefail
shopt -s nullglob

log "Starting RAT (read-level classification)"

# ------------------------------------------------------------
# RAT is extremely memory hungry
# ------------------------------------------------------------
THREADS="${THREADS:-8}"

export OMP_NUM_THREADS="$THREADS"

TRIM_DIR="$OUTPUT_DIR/trimmed"

# TIARA contigs
TIARA_DIR="$OUTPUT_DIR/tiara"

# Cleaned contigs
TIARA_CLEAN_DIR="$OUTPUT_DIR/tiara_clean"

# Existing CAT outputs
CAT_CONTIG_DIR="$OUTPUT_DIR/cat_contigs"

RAT_DIR="$OUTPUT_DIR/rat_reads"
LOG_DIR="$OUTPUT_DIR/logs/rat_reads"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$RAT_DIR" "$LOG_DIR" "$TIARA_CLEAN_DIR"

# ============================================================
# CHECK DEPENDENCIES
# ============================================================
command -v "$CAT_CMD" >/dev/null || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null || {
    log "ERROR: diamond not found"
    exit 1
}

command -v bwa >/dev/null || {
    log "ERROR: bwa not found"
    exit 1
}

command -v samtools >/dev/null || {
    log "ERROR: samtools not found"
    exit 1
}

# ============================================================
# COLLECT SAMPLES TO RUN
# ============================================================
PENDING_SAMPLES=()

for f in "$TRIM_DIR"/*_trimmed_R1.fastq.gz; do

    [[ -s "$f" ]] || continue

    SAMPLE=$(basename "$f" | sed 's/_trimmed_R1.fastq.gz//')

    # --------------------------------------------------------
    # OPTIONAL: explicit target sample
    # --------------------------------------------------------
    if [[ -n "${TARGET_SAMPLE:-}" ]]; then

        if [[ "$SAMPLE" != "$TARGET_SAMPLE" ]]; then
            continue
        fi
    fi

    PREFIX="$RAT_DIR/$SAMPLE/${SAMPLE}"

    # --------------------------------------------------------
    # SUCCESS CHECK
    # --------------------------------------------------------
    COMPLETE_FILE="${PREFIX}.RAT.complete.abundance.txt"
    CONTIG_FILE="${PREFIX}.RAT.contig.abundance.txt"
    READ_FILE="${PREFIX}.RAT.read2classification.txt"

    if [[ -s "$COMPLETE_FILE" && -s "$CONTIG_FILE" && -s "$READ_FILE" ]]; then

        if [[ "${FORCE:-0}" != "1" ]]; then
            log "Already completed → skipping $SAMPLE"
            continue
        fi
    fi

    PENDING_SAMPLES+=("$SAMPLE")
done

N=${#PENDING_SAMPLES[@]}

log "Samples pending RAT: $N"

[[ "$N" -gt 0 ]] || {
    log "No samples need RAT"
    exit 0
}

printf '%s\n' "${PENDING_SAMPLES[@]}" \
    > "$LOG_DIR/rat_pending_samples.txt"

log "Pending samples:"
cat "$LOG_DIR/rat_pending_samples.txt"

# ============================================================
# SAMPLE SELECTION
# ============================================================

# ------------------------------------------------------------
# CASE 1: explicit TARGET_SAMPLE
# ------------------------------------------------------------
if [[ -n "${TARGET_SAMPLE:-}" ]]; then

    SAMPLE="$TARGET_SAMPLE"

# ------------------------------------------------------------
# CASE 2: SLURM ARRAY
# ------------------------------------------------------------
else

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        log "ERROR: Must run as SLURM array or use TARGET_SAMPLE"
        exit 1
    fi

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then
        log "No sample for this array index"
        exit 0
    fi

    SAMPLE="${PENDING_SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
fi

log "Selected sample: $SAMPLE"

# ============================================================
# INPUT FILES
# ============================================================
R1="$TRIM_DIR/${SAMPLE}_trimmed_R1.fastq.gz"
R2="$TRIM_DIR/${SAMPLE}_trimmed_R2.fastq.gz"

RAW_CONTIGS="$TIARA_DIR/$SAMPLE/euk.fa"

CLEAN_DIR="$TIARA_CLEAN_DIR/$SAMPLE"
mkdir -p "$CLEAN_DIR"

CONTIGS="$CLEAN_DIR/euk_clean.fa"

RAW_C2C="$CAT_CONTIG_DIR/$SAMPLE/${SAMPLE}.contig2classification.txt"

CAT_C2C="$CAT_CONTIG_DIR/$SAMPLE/${SAMPLE}.contig2classification.clean.txt"

OUTDIR="$RAT_DIR/$SAMPLE"
PREFIX="$OUTDIR/${SAMPLE}"

TMPDIR="$OUTDIR/tmp"

DUMMY_U2C="$OUTDIR/dummy_u2c.txt"

log "----------------------------------------"
log "Sample:         $SAMPLE"
log "R1:             $R1"
log "R2:             $R2"
log "Raw contigs:    $RAW_CONTIGS"
log "Clean contigs:  $CONTIGS"
log "CAT C2C:        $CAT_C2C"
log "TMPDIR:         $TMPDIR"
log "THREADS:        $THREADS"

# ============================================================
# INPUT CHECKS
# ============================================================
[[ -s "$R1" ]] || {
    log "ERROR: Missing R1"
    exit 1
}

[[ -s "$R2" ]] || {
    log "ERROR: Missing R2"
    exit 1
}

[[ -s "$RAW_CONTIGS" ]] || {
    log "ERROR: Missing Tiara euk.fa"
    exit 1
}

[[ -s "$RAW_C2C" ]] || {
    log "ERROR: Missing CAT contig2classification"
    exit 1
}

[[ -d "$CAT_TAX" ]] || {
    log "ERROR: Missing CAT taxonomy"
    exit 1
}

# ============================================================
# CLEAN OUTPUT ONLY FOR CURRENT SAMPLE
# ============================================================
rm -rf "$OUTDIR"

mkdir -p "$OUTDIR" "$TMPDIR"

touch "$DUMMY_U2C"

# ============================================================
# CLEAN FASTA HEADERS
# ============================================================
log "Cleaning FASTA headers"

awk '
/^>/ {
    sub(/^>/,"")
    split($0,a," ")
    print ">"a[1]
    next
}
{print}
' \
"$RAW_CONTIGS" \
> "$CONTIGS"

# ============================================================
# CLEAN CAT contig2classification
# ============================================================
log "Cleaning CAT contig2classification"

awk '
BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
    split($1,a," ")
    $1=a[1]
    print
}
' \
"$RAW_C2C" \
> "$CAT_C2C"

# ============================================================
# VALIDATE CLEAN FILES
# ============================================================
[[ -s "$CONTIGS" ]] || {
    log "ERROR: Cleaned contigs file failed"
    exit 1
}

[[ -s "$CAT_C2C" ]] || {
    log "ERROR: Cleaned CAT file failed"
    exit 1
}

# ============================================================
# RUN RAT
# ============================================================
log "Running RAT"

"$CAT_CMD" reads \
    --mode cr \
    -c "$CONTIGS" \
    --c2c "$CAT_C2C" \
    --u2c "$DUMMY_U2C" \
    -1 "$R1" \
    -2 "$R2" \
    -t "$CAT_TAX" \
    -o "$PREFIX" \
    --tmpdir "$TMPDIR" \
    --read2classification \
    --force \
    --nproc "$THREADS" \
    2>&1 | tee "$LOG_DIR/${SAMPLE}.log"

# ============================================================
# VALIDATE OUTPUT
# ============================================================
[[ -s "${PREFIX}.RAT.read2classification.txt" ]] || {
    log "ERROR: Missing read2classification output"
    exit 1
}

[[ -s "${PREFIX}.RAT.complete.abundance.txt" ]] || {
    log "ERROR: Missing complete abundance output"
    exit 1
}

[[ -s "${PREFIX}.RAT.contig.abundance.txt" ]] || {
    log "ERROR: Missing contig abundance output"
    exit 1
}

# ============================================================
# OUTPUT SUMMARY
# ============================================================
log "Generated outputs:"
ls -lh "$OUTDIR"/*.RAT.* || true

# ============================================================
# CLEANUP
# ============================================================
rm -f "$DUMMY_U2C"

log "Finished RAT for $SAMPLE"
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
### RAT (READ-LEVEL — TIARA PROK CONTIGS + CAT PROK OUTPUTS)
### RUN NEW SAMPLES + RERUN FAILED ONLY
# ============================================================
rat_reads_prok)

set -euo pipefail
shopt -s nullglob

log "Starting RAT (read-level classification, prok)"

# ============================================================
# CONFIG
# ============================================================
THREADS="${THREADS:-8}"

export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

TRIM_DIR="$OUTPUT_DIR/trimmed"

# Tiara outputs
TIARA_DIR="$OUTPUT_DIR/tiara"

# Cleaned prok contigs
TIARA_CLEAN_DIR="$OUTPUT_DIR/tiara_clean_prok"

# Existing CAT prok outputs
CAT_CONTIG_DIR="$OUTPUT_DIR/cat_contigs_prok"

RAT_DIR="$OUTPUT_DIR/rat_reads_prok"

LOG_DIR="$OUTPUT_DIR/logs/rat_reads_prok"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"

CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"

CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$RAT_DIR" "$LOG_DIR" "$TIARA_CLEAN_DIR"

# ============================================================
# CHECK DEPENDENCIES
# ============================================================
command -v "$CAT_CMD" >/dev/null || {
    log "ERROR: CAT_pack not found"
    exit 1
}

command -v diamond >/dev/null || {
    log "ERROR: diamond not found"
    exit 1
}

command -v bwa >/dev/null || {
    log "ERROR: bwa not found"
    exit 1
}

command -v samtools >/dev/null || {
    log "ERROR: samtools not found"
    exit 1
}

# ============================================================
# COLLECT PENDING SAMPLES
# ============================================================
PENDING_SAMPLES=()

for f in "$TRIM_DIR"/*_trimmed_R1.fastq.gz; do

    [[ -s "$f" ]] || continue

    SAMPLE=$(basename "$f" | sed 's/_trimmed_R1.fastq.gz//')

    # --------------------------------------------------------
    # OPTIONAL TARGET SAMPLE
    # --------------------------------------------------------
    if [[ -n "${TARGET_SAMPLE:-}" ]]; then

        [[ "$SAMPLE" == "$TARGET_SAMPLE" ]] || continue
    fi

    PREFIX="$RAT_DIR/$SAMPLE/${SAMPLE}"

    COMPLETE_FILE="${PREFIX}.RAT.complete.abundance.txt"

    CONTIG_FILE="${PREFIX}.RAT.contig.abundance.txt"

    READ_FILE="${PREFIX}.RAT.read2classification.txt"

    # --------------------------------------------------------
    # SKIP COMPLETED
    # --------------------------------------------------------
    if [[ -s "$COMPLETE_FILE" &&
          -s "$CONTIG_FILE" &&
          -s "$READ_FILE" &&
          "${FORCE:-0}" != "1" ]]; then

        log "Already completed → skipping $SAMPLE"

        continue
    fi

    PENDING_SAMPLES+=("$SAMPLE")
done

N=${#PENDING_SAMPLES[@]}

log "Samples pending RAT: $N"

[[ "$N" -gt 0 ]] || {

    log "No samples require RAT"

    exit 0
}

printf '%s\n' "${PENDING_SAMPLES[@]}" \
    > "$LOG_DIR/rat_prok_pending_samples.txt"

log "Pending samples:"
cat "$LOG_DIR/rat_prok_pending_samples.txt"

# ============================================================
# SAMPLE SELECTION
# ============================================================

# ------------------------------------------------------------
# TARGET SAMPLE
# ------------------------------------------------------------
if [[ -n "${TARGET_SAMPLE:-}" ]]; then

    SAMPLE="$TARGET_SAMPLE"

# ------------------------------------------------------------
# ARRAY MODE
# ------------------------------------------------------------
else

    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {

        log "ERROR: Must run as array job or TARGET_SAMPLE"

        exit 1
    }

    if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > N )); then

        log "No sample assigned to this array index"

        exit 0
    fi

    SAMPLE="${PENDING_SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
fi

log "Selected sample: $SAMPLE"

# ============================================================
# INPUT FILES
# ============================================================
R1="$TRIM_DIR/${SAMPLE}_trimmed_R1.fastq.gz"

R2="$TRIM_DIR/${SAMPLE}_trimmed_R2.fastq.gz"

RAW_CONTIGS="$TIARA_DIR/$SAMPLE/prok.fa"

CLEAN_DIR="$TIARA_CLEAN_DIR/$SAMPLE"

mkdir -p "$CLEAN_DIR"

CONTIGS="$CLEAN_DIR/prok_clean.fa"

RAW_C2C="$CAT_CONTIG_DIR/$SAMPLE/${SAMPLE}.contig2classification.txt"

CAT_C2C="$CAT_CONTIG_DIR/$SAMPLE/${SAMPLE}.contig2classification.clean.txt"

OUTDIR="$RAT_DIR/$SAMPLE"

PREFIX="$OUTDIR/${SAMPLE}"

TMPDIR="$OUTDIR/tmp"

DUMMY_U2C="$OUTDIR/dummy_u2c.txt"

# ============================================================
# DEBUG
# ============================================================
log "----------------------------------------"

log "Sample:         $SAMPLE"

log "R1:             $R1"

log "R2:             $R2"

log "Raw contigs:    $RAW_CONTIGS"

log "Clean contigs:  $CONTIGS"

log "CAT C2C:        $CAT_C2C"

log "TMPDIR:         $TMPDIR"

log "THREADS:        $THREADS"

# ============================================================
# INPUT CHECKS
# ============================================================
[[ -s "$R1" ]] || {
    log "ERROR: Missing R1"
    exit 1
}

[[ -s "$R2" ]] || {
    log "ERROR: Missing R2"
    exit 1
}

[[ -s "$RAW_CONTIGS" ]] || {
    log "ERROR: Missing Tiara prok.fa"
    exit 1
}

[[ -s "$RAW_C2C" ]] || {
    log "ERROR: Missing CAT contig2classification"
    exit 1
}

[[ -d "$CAT_TAX" ]] || {
    log "ERROR: Missing CAT taxonomy"
    exit 1
}

# ============================================================
# CLEAN CURRENT OUTPUT
# ============================================================
rm -rf "$OUTDIR"

mkdir -p "$OUTDIR" "$TMPDIR"

touch "$DUMMY_U2C"

# ============================================================
# CLEAN FASTA HEADERS
# ============================================================
log "Cleaning FASTA headers"

awk '
/^>/ {
    sub(/^>/,"")
    split($0,a," ")
    print ">"a[1]
    next
}
{print}
' \
"$RAW_CONTIGS" \
> "$CONTIGS"

# ============================================================
# CLEAN contig2classification
# ============================================================
log "Cleaning contig2classification"

awk '
BEGIN{FS=OFS="\t"}

/^#/ {
    print
    next
}

{
    split($1,a," ")

    $1=a[1]

    print
}
' \
"$RAW_C2C" \
> "$CAT_C2C"

# ============================================================
# VALIDATE CLEAN FILES
# ============================================================
[[ -s "$CONTIGS" ]] || {
    log "ERROR: Failed to create cleaned prok contigs"
    exit 1
}

[[ -s "$CAT_C2C" ]] || {
    log "ERROR: Failed to create cleaned contig2classification"
    exit 1
}

# ============================================================
# RUN RAT
# ============================================================
log "Running RAT (prok)"

"$CAT_CMD" reads \
    --mode cr \
    -c "$CONTIGS" \
    --c2c "$CAT_C2C" \
    --u2c "$DUMMY_U2C" \
    -1 "$R1" \
    -2 "$R2" \
    -t "$CAT_TAX" \
    -o "$PREFIX" \
    --tmpdir "$TMPDIR" \
    --read2classification \
    --force \
    --nproc "$THREADS" \
    2>&1 | tee "$LOG_DIR/${SAMPLE}.log"

# ============================================================
# VALIDATE OUTPUT
# ============================================================
[[ -s "${PREFIX}.RAT.read2classification.txt" ]] || {
    log "ERROR: Missing read2classification output"
    exit 1
}

[[ -s "${PREFIX}.RAT.complete.abundance.txt" ]] || {
    log "ERROR: Missing complete abundance output"
    exit 1
}

[[ -s "${PREFIX}.RAT.contig.abundance.txt" ]] || {
    log "ERROR: Missing contig abundance output"
    exit 1
}

# ============================================================
# OUTPUT SUMMARY
# ============================================================
log "Generated outputs:"

ls -lh "$OUTDIR"/*.RAT.* || true

# ============================================================
# CLEANUP
# ============================================================
rm -f "$DUMMY_U2C"

log "Finished RAT for $SAMPLE"

;;

# ============================================================
### CAT (CONTIG-LEVEL — TIARA EUKARYOTIC CONTIGS)
# ============================================================
cat_contigs)

set -euo pipefail
shopt -s nullglob

log "Starting CAT (contig-level, euk)"

THREADS="${THREADS:-16}"
TARGET_SAMPLE="${TARGET_SAMPLE:-}"

TIARA_DIR="$OUTPUT_DIR/tiara"
CAT_DIR="$OUTPUT_DIR/cat_contigs"
LOG_DIR="$OUTPUT_DIR/logs/cat_contigs"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$CAT_DIR" "$LOG_DIR"

# -------------------- DEPENDENCY CHECK --------------------
command -v "$CAT_CMD" >/dev/null || { log "ERROR: CAT_pack not found"; exit 1; }
command -v diamond >/dev/null || { log "ERROR: DIAMOND not found"; exit 1; }

# ============================================================
# FUNCTION TO PROCESS ONE SAMPLE
# ============================================================
run_cat() {

    local SAMPLE="$1"
    local COND="$TIARA_DIR/$SAMPLE"

    CONTIGS="$COND/euk.fa"
    OUTDIR="$CAT_DIR/$SAMPLE"
    PREFIX="$OUTDIR/${SAMPLE}"

    log "----------------------------------------"
    log "Processing: $SAMPLE"

    [[ -s "$CONTIGS" ]] || {
        log "No euk.fa → skip"
        return
    }

    # -------------------- SKIP IF DONE --------------------
    if [[ -s "${PREFIX}.contigs2classification.txt" && "${FORCE:-0}" != "1" ]]; then
        log "Already done → skipping"
        return
    fi

    # -------------------- CLEAN --------------------
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"

    # -------------------- RUN CAT --------------------
    log "Running CAT_pack"

    "$CAT_CMD" contigs \
        -c "$CONTIGS" \
        -d "$CAT_DB" \
        -t "$CAT_TAX" \
        -o "$PREFIX" \
        --nproc "$THREADS" \
        2>&1 | tee "$LOG_DIR/${SAMPLE}.log"

    [[ -s "${PREFIX}.contigs2classification.txt" ]] || {
        log "ERROR: CAT failed"
        exit 1
    }

    log "Finished: $SAMPLE"
}

# ============================================================
# MODE SELECTION
# ============================================================

# -------------------- MANUAL MODE --------------------
if [[ -n "$TARGET_SAMPLE" ]]; then
    log "Manual mode → $TARGET_SAMPLE"
    run_cat "$TARGET_SAMPLE"

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
    run_cat "$SAMPLE"

# -------------------- AUTO MODE (RUN ONLY MISSING) --------------------
else
    log "Auto mode → running only missing samples"

    for d in "$TIARA_DIR"/*; do
        [[ -d "$d" ]] || continue
        SAMPLE=$(basename "$d")

        CONTIGS="$d/euk.fa"
        PREFIX="$CAT_DIR/$SAMPLE/${SAMPLE}"

        if [[ -s "$CONTIGS" ]] && [[ ! -s "${PREFIX}.contigs2classification.txt" ]]; then
            run_cat "$SAMPLE"
        else
            log "Skipping $SAMPLE (done or no data)"
        fi
    done
fi

log "CAT step completed"
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
# ============================================================
cat_mito_contigs)

set -euo pipefail
shopt -s nullglob

log "Starting CAT (mitochondrial contigs) step"

THREADS="${THREADS:-16}"

TIARA_DIR="$OUTPUT_DIR/tiara"
CAT_DIR="$OUTPUT_DIR/cat_mito"
LOG_DIR="$OUTPUT_DIR/logs/cat_mito"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$CAT_DIR" "$LOG_DIR"

# -------------------- DEPENDENCY CHECK --------------------
command -v "$CAT_CMD" >/dev/null || { log "ERROR: CAT_pack not found"; exit 1; }
command -v diamond >/dev/null || { log "ERROR: DIAMOND not found"; exit 1; }

# -------------------- COLLECT SAMPLES --------------------
mapfile -t SAMPLES < <(
    for d in "$TIARA_DIR"/*; do
        [[ -d "$d" ]] || continue
        [[ -s "$d/mito.fa" ]] && echo "$d"
    done | sort
)

N=${#SAMPLES[@]}
log "Samples with mito contigs: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No mito.fa found"; exit 1; }

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID > N )); then
    log "No sample for this array index → exiting"
    exit 0
fi

# -------------------- SELECT SAMPLE --------------------
COND="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$COND")

CONTIGS="$COND/mito.fa"
OUTDIR="$CAT_DIR/$SAMPLE"
PREFIX="$OUTDIR/${SAMPLE}_mito"

log "----------------------------------------"
log "Processing: $SAMPLE"
log "Input (mito): $CONTIGS"

# -------------------- SKIP --------------------
if [[ -s "${PREFIX}.contigs2classification.txt" && "${FORCE:-0}" != "1" ]]; then
    log "CAT already done → skipping"
    exit 0
fi

# -------------------- CLEAN --------------------
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# -------------------- RUN CAT --------------------
log "Running CAT_pack (mito contigs)"

"$CAT_CMD" contigs \
    -c "$CONTIGS" \
    -d "$CAT_DB" \
    -t "$CAT_TAX" \
    -o "$PREFIX" \
    --nproc "$THREADS" \
    2>&1 | tee "$LOG_DIR/${SAMPLE}_mito.log"

# -------------------- VALIDATE --------------------
[[ -s "${PREFIX}.contigs2classification.txt" ]] || {
    log "ERROR: CAT failed"
    exit 1
}

log "Finished CAT (mito) for $SAMPLE"
;;

# ============================================================
### BAT (BIN-LEVEL — REMAG BINS)
# ============================================================

bat_bins)

set -euo pipefail
shopt -s nullglob

log "Starting BAT (bin-level) step"

THREADS="${THREADS:-16}"

REMAG_DIR="$OUTPUT_DIR/remag"
BAT_DIR="$OUTPUT_DIR/bat_bins"
LOG_DIR="$OUTPUT_DIR/logs/bat_bins"

CAT_DB="${CAT_DB:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/db}"
CAT_TAX="${CAT_TAX:-/staging/leuven/stg_00100/Rodgee/databases/CAT/20241212_CAT_nr_website/tax}"
CAT_CMD="${CAT_CMD:-CAT_pack}"

mkdir -p "$BAT_DIR" "$LOG_DIR"

# -------------------- DEPENDENCY CHECK --------------------
command -v "$CAT_CMD" >/dev/null || { log "ERROR: CAT_pack not found"; exit 1; }
command -v diamond >/dev/null || { log "ERROR: DIAMOND not found"; exit 1; }

# -------------------- COLLECT SAMPLES --------------------
mapfile -t SAMPLES < <(
    for d in "$REMAG_DIR"/*; do
        [[ -d "$d" ]] && echo "$d"
    done | sort
)

N=${#SAMPLES[@]}
log "Total samples: $N"

[[ "$N" -gt 0 ]] || { log "ERROR: No reMAG samples found"; exit 1; }

# -------------------- ARRAY CHECK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log "ERROR: Run as array job"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID > N )); then
    log "No sample for this array index → exiting"
    exit 0
fi

# -------------------- SELECT SAMPLE --------------------
SAMPLE_DIR="${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
SAMPLE=$(basename "$SAMPLE_DIR")

BINS_DIR="$SAMPLE_DIR/bins"
OUTDIR="$BAT_DIR/$SAMPLE"
PREFIX="$OUTDIR/${SAMPLE}"

log "----------------------------------------"
log "Processing: $SAMPLE"
log "Bins dir: $BINS_DIR"

# -------------------- CHECK INPUT --------------------
if [[ ! -d "$BINS_DIR" ]] || ! compgen -G "$BINS_DIR/*.fa" > /dev/null; then
    log "No bins found → skipping"
    exit 0
fi

# -------------------- SKIP --------------------
if [[ -s "${PREFIX}.bin2classification.txt" && "${FORCE:-0}" != "1" ]]; then
    log "BAT already done → skipping"
    exit 0
fi

# -------------------- CLEAN --------------------
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# -------------------- RUN BAT --------------------
log "Running CAT_pack (bins)"
log "Output prefix: $PREFIX"

"$CAT_CMD" bins \
    -b "$BINS_DIR" \
    -d "$CAT_DB" \
    -t "$CAT_TAX" \
    -o "$PREFIX" \
    --bin_suffix .fa \
    --nproc "$THREADS" \
    2>&1 | tee "$LOG_DIR/${SAMPLE}.log"

# -------------------- VALIDATE --------------------
[[ -s "${PREFIX}.bin2classification.txt" ]] || {
    log "ERROR: BAT failed"
    exit 1
}

log "Finished BAT for $SAMPLE"
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

log "Creating merged summary table"

FIRST=1

for FILE in "$SUMMARY_DIR"/*.taxonomy_structured.tsv; do

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
        ' "$FILE" > "$MERGED"

        FIRST=0

    else

        awk -v sample="$SAMPLE" '
        BEGIN { OFS="\t" }

        NR>1 {
            print sample,$0
        }
        ' "$FILE" >> "$MERGED"
    fi

done

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
    COV_FILE="$OUTPUT_DIR/coverage/${SAMPLE}.coverage.tsv"

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
# COMPARE INCUBATION BASED ON BINS: TAXONOMY BINS + BUSCO
# ============================================================
compare_incubation_full_taxonomy)

set -euo pipefail
shopt -s nullglob

log() { echo "[$(date '+%F %T')] $*"; }

SUMMARY_DIR="$OUTPUT_DIR/bat_bins_summary"
META="/staging/leuven/stg_00100/Rodgee/caves_data/metadata/biological_metadata.tsv"
BUSCO="$OUTPUT_DIR/busco_classified/busco_summary.tsv"

OUT_BASE="$OUTPUT_DIR/comparison_incubation/full_taxonomy.tsv"
OUT_FINAL="$OUTPUT_DIR/comparison_incubation/full_taxonomy_busco.tsv"

mkdir -p "$(dirname "$OUT_BASE")"

# STEP 1: BUILD TAXONOMY TABLE (ONLY IF MISSING)
if [[ -s "$OUT_BASE" && "${FORCE:-0}" != "1" ]]; then
    log "Full taxonomy table exists → skipping generation"
else
    log "Generating full taxonomy table"

    awk -v meta="$META" '
    BEGIN {
        FS="\t"; OFS="\t"

        while ((getline < meta) > 0) {
            if ($1=="GC_code") continue
            incubation[$1]=$2
            samplecode[$1]=$3
        }
        close(meta)

        print "sample","sample_code","incubation","bin","superkingdom","kingdom","phylum","class","order","family","genus","species"
    }

    FNR==1 { next }

    {
        n = split(FILENAME, a, "/")
        file = a[n]
        split(file, b, "\\.")
        sample = b[1]

        inc  = incubation[sample]
        code = samplecode[sample]

        if (inc == "") next

        print sample, code, inc, $1, $2, $3, $4, $5, $6, $7, $8, $9
    }
    ' "$SUMMARY_DIR"/S*.taxonomy_structured.tsv > "$OUT_BASE"

    log "Created: $OUT_BASE"
fi

# STEP 2: MERGE BUSCO (STRICT sample|bin MATCH)
if [[ ! -s "$BUSCO" ]]; then
    log "WARNING: BUSCO summary not found → skipping merge"
    exit 0
fi

if [[ -s "$OUT_FINAL" && "${FORCE:-0}" != "1" ]]; then
    log "Final taxonomy+BUSCO exists → skipping merge"
    exit 0
fi

log "Preparing composite keys for join"

# -------------------- PREP TAX --------------------
awk '
NR>1 {
    sample = $1
    bin = $4

    # Fix: remove .fa
    sub(/\.fa$/, "", bin)

    key = sample "|" bin

    print key"\t"$0
}
' "$OUT_BASE" \
| sort -k1,1 > tmp.tax

# -------------------- PREP BUSCO --------------------
awk '
NR>1 {
    sample = $1
    bin = $2

    key = sample "|" bin

    print key"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8
}
' "$BUSCO" \
| sort -k1,1 > tmp.busco

# JOIN

join -t $'\t' tmp.tax tmp.busco \
| awk '
BEGIN {
    OFS="\t"
    print "sample","sample_code","incubation","bin","superkingdom","kingdom","phylum","class","order","family","genus","species","completeness","single","duplicated","fragmented","missing","category"
}
{
    # Layout after join:
    # $1 = key
    # $2-$13 = taxonomy columns
    # $14-$19 = BUSCO columns

    sample       = $2
    sample_code  = $3
    incubation   = $4
    bin          = $5

    print sample, sample_code, incubation, bin, \
          $6,$7,$8,$9,$10,$11,$12,$13, \
          $14,$15,$16,$17,$18,$19
}
' > "$OUT_FINAL"

# -------------------- CLEAN --------------------
rm -f tmp.tax tmp.busco

log "Created: $OUT_FINAL"
;;

# ============================================================
# COMPARISON OF INCUBATION BASED ON BINS: QUAST + BUSCO
# ============================================================
compare_incubation_quast_busco)

set -euo pipefail

log() { echo "[$(date '+%F %T')] $*"; }

IN="$OUTPUT_DIR/quast/quast_busco_merged.tsv"
OUT="$OUTPUT_DIR/quast/quast_busco_comparison.tsv"

log "Creating combined comparison (final clean)"

{
awk '
BEGIN {
    FS="\t"; OFS="\t"
}

NR==1 { next }

{
    base = $2
    sub(/-[LS]$/, "", base)

    if ($3=="long") long[base]=$0
    else if ($3=="short") short[base]=$0
}

END {
    for (k in long) {
        if (k in short) {

            split(long[k], L, "\t")
            split(short[k], S, "\t")

            # numeric conversion
            n50_L = L[6]+0; n50_S = S[6]+0
            len_L = L[5]+0; len_S = S[5]+0
            comp_L = L[8]+0; comp_S = S[8]+0
            frag_L = L[9]+0; frag_S = S[9]+0
            gc_L = L[7]+0; gc_S = S[7]+0
            cont_L = L[4]+0; cont_S = S[4]+0
            bins_L = L[11]+0; bins_S = S[11]+0

            # FIX: use IF instead of next
            if (!(bins_L == 0 && bins_S == 0)) {

                printf "%s\t%s\t%s\t%d\t%d\t%.3f\t%.3f\t%.2f\t%d\t%d\t%d\n",
                       k,
                       L[1], S[1],
                       (n50_L - n50_S),
                       (len_L - len_S),
                       (comp_L - comp_S),
                       (frag_L - frag_S),
                       (gc_L - gc_S),
                       (cont_L - cont_S),
                       bins_L, bins_S
            }
        }
    }
}
' "$IN" | sort
} | awk '
BEGIN {
    OFS="\t"
    print "sample_code",
          "long_sample","short_sample",
          "delta_N50","delta_length",
          "delta_completeness","delta_fragmented",
          "delta_GC","delta_contigs",
          "long_bins","short_bins"
}
{ print }
' > "$OUT"

log "Created: $OUT"

;;

# ============================================================
# COMPARISON BASED ON BINS : QUAST + BUSCO + COVERAGE
# ============================================================
compare_incubation_full)

set -euo pipefail

log() { echo "[$(date '+%F %T')] $*"; }

QB="$OUTPUT_DIR/quast/quast_busco_merged.tsv"
COV="$OUTPUT_DIR/coverage/coverage_per_sample.tsv"
OUT="$OUTPUT_DIR/quast/full_comparison_quast_busco_cov.tsv"

log "Creating full comparison (FINAL CLEAN)"

# -------------------- PREP --------------------
tail -n +2 "$QB"  | sort -k1,1 > tmp.qb
tail -n +2 "$COV" | sort -k1,1 > tmp.cov

# explicit join (correct columns)
join -t $'\t' \
-o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.2,2.3,2.4 \
tmp.qb tmp.cov > tmp.all

# -------------------- BUILD TABLE --------------------
awk '
BEGIN {
    FS="\t"; OFS="\t"
}

{
    code = $2
    incub = $3

    base = code
    sub(/-[LS]$/, "", base)

    if (incub=="long")  long[base]=$0
    if (incub=="short") short[base]=$0
}

END {

    for (k in long) {
        if (k in short) {

            split(long[k], L, "\t")
            split(short[k], S, "\t")

            # parse
            n50_L=L[6]+0; n50_S=S[6]+0
            len_L=L[5]+0; len_S=S[5]+0
            comp_L=L[8]+0; comp_S=S[8]+0
            frag_L=L[9]+0; frag_S=S[9]+0
            bins_L=L[11]+0; bins_S=S[11]+0
            cov_tot_L=L[12]+0; cov_tot_S=S[12]+0
            cov_mean_L=L[14]+0; cov_mean_S=S[14]+0

            if (bins_L==0 && bins_S==0) continue

            printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%.3f\t%.3f\t%.0f\t%.0f\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.0f\t%d\n",
                   k,
                   L[1], S[1],
                   n50_L, n50_S,
                   len_L, len_S,
                   comp_L, comp_S,
                   frag_L, frag_S,
                   bins_L, bins_S,
                   cov_mean_L, cov_mean_S,
                   cov_tot_L, cov_tot_S,
                   (n50_L-n50_S),
                   (len_L-len_S),
                   (comp_L-comp_S),
                   (frag_L-frag_S),
                   (cov_mean_L-cov_mean_S),
                   (cov_tot_L-cov_tot_S),
                   (bins_L-bins_S)
        }
    }
}
' tmp.all | sort > tmp.body

# -------------------- ADD HEADER (AFTER SORT) --------------------
{
echo -e "sample_code\tlong_sample\tshort_sample\tlong_N50\tshort_N50\tlong_length\tshort_length\tlong_completeness\tshort_completeness\tlong_fragmented\tshort_fragmented\tlong_bins\tshort_bins\tlong_mean_cov\tshort_mean_cov\tlong_total_cov\tshort_total_cov\tdelta_N50\tdelta_length\tdelta_completeness\tdelta_fragmented\tdelta_mean_cov\tdelta_total_cov\tdelta_bins"
cat tmp.body
} > "$OUT"

rm -f tmp.qb tmp.cov tmp.all tmp.body

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