#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=0:10:00
#$ -l h_vmem=2G
#$ -N collect_ranked_models

# === CONFIGURATION ===
PROT="$1"  # e.g. trimmed_lysozyme
if [[ -z "$PROT" ]]; then
    echo "[!] Usage: qsub collect_ranked_models.sh <PROT>"
    exit 1
fi

SOURCE_BASE="/data/home/ap24012/scratch_preserved/af2_results"
DEST_DIR="/data/home/ap24012/scratch_preserved/af2_evaluation/predictions/${PROT}"

mkdir -p "$DEST_DIR"

echo "[*] Collecting ranked models for: $PROT"
echo "[*] Output directory: $DEST_DIR"

# === Main Loop ===
find "$SOURCE_BASE" -type f -name "ranked_*.pdb" | grep "$PROT" | while read -r file; do
    CLUSTER=$(echo "$file" | grep -o "${PROT}_U[0-9]\{2,3\}-[0-9]\{3\}")
    BASENAME=$(basename "$file")
    NEWNAME="${CLUSTER}_${BASENAME}"
    cp "$file" "${DEST_DIR}/${NEWNAME}"
    echo "  [+] Copied $file → ${NEWNAME}"
done

echo "[✓] All ranked models saved to: $DEST_DIR"
