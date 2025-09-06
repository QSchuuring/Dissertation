#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=1:00:00
#$ -l h_vmem=8G
#$ -N cluster_pipeline

# === CONFIGURATION ===
PROT="$1"
BASE_SCRATCH="/data/home/ap24012/scratch_preserved/af2_runs/${PROT}"
RAW_STO_DIR="/data/home/ap24012/alphafold_input/${PROT}/msas_raw_sto"
SCRIPT_DIR="/data/home/ap24012/scripts"
MERGED_A3M_DIR="${BASE_SCRATCH}/merged_a3m"
CLUSTER_OUT="${BASE_SCRATCH}/clusters"
ENV_PATH="/data/home/ap24012/.conda/envs/afcluster-env"
PERL_PATH="/share/apps/rocky9/spack/apps/linux-rocky9-x86_64_v4/gcc-12.2.0/perl/5.38.0-hkeva47/bin/perl"

mkdir -p "$MERGED_A3M_DIR" "$CLUSTER_OUT"

# === STEP 0: ENVIRONMENT ===
module load miniforge perl mafft
source /share/apps/rocky9/general/apps/miniforge/24.7.1/etc/profile.d/conda.sh
conda activate $ENV_PATH
export PATH="$CONDA_PREFIX/bin:$PATH"

echo "[*] Conda: $CONDA_PREFIX"
echo "[*] Python: $(which python)"
echo "[*] Perl: $PERL_PATH"

# === STEP 1: Convert .sto → .a3m (and fix to single-line format) ===
echo "[*] Converting raw .sto to single-line .a3m..."
for file in "$RAW_STO_DIR"/*.sto; do
    base=$(basename "$file" .sto)
    out_a3m="${MERGED_A3M_DIR}/${base}.a3m"

    $PERL_PATH "${SCRIPT_DIR}/reformat.pl" sto a3m "$file" "$out_a3m"

    # Fix multiline sequences to single-line A3M
    awk '/^>/{if (seq) print seq; print; seq=""; next} {seq=seq $0} END{print seq}' \
      "$out_a3m" > "${out_a3m}.tmp" && mv "${out_a3m}.tmp" "$out_a3m"
done
echo "[✓] All .sto converted to fixed .a3m"

# === STEP 2: Merge all A3Ms ===
MERGED_A3M="${MERGED_A3M_DIR}/${PROT}.a3m"
rm -f "$MERGED_A3M"
cat "${MERGED_A3M_DIR}/"*.a3m > "$MERGED_A3M"
echo "[✓] Merged A3Ms → $MERGED_A3M"

# === STEP 3: Run AFCluster ===
echo "[*] Running clustering..."
python /data/home/ap24012/programs/AF_Cluster/scripts/ClusterMSA.py \
  "$PROT" \
  -i "$MERGED_A3M" \
  -o "$CLUSTER_OUT" \
  --run_PCA
echo "[✓] Clustering complete → $CLUSTER_OUT"

# === STEP 4: For each cluster A3M, realign and create .sto files ===
echo "[*] Processing clusters for AF2..."

for cluster_a3m in ${CLUSTER_OUT}/${PROT}_U{10..30..2}-*.a3m; do
    CLUSTER_NAME=$(basename "$cluster_a3m" .a3m)
    OUTDIR="${BASE_SCRATCH}/${CLUSTER_NAME}"
    mkdir -p "$OUTDIR"

    echo "  - Aligning $CLUSTER_NAME"

    # Realign with MAFFT
    mafft --auto "$cluster_a3m" > "${OUTDIR}/${CLUSTER_NAME}_aligned.a3m"

    # Postprocess
    python3 "${SCRIPT_DIR}/postprocess_a3m.py" \
        "${OUTDIR}/${CLUSTER_NAME}_aligned.a3m" \
        "${OUTDIR}/${CLUSTER_NAME}_aligned.a3m"

    # Reformat to .sto and fix (for each of the 3 databases)
    for db in uniref90_hits small_bfd_hits mgnify_hits; do
        a3m_path="${OUTDIR}/${CLUSTER_NAME}_aligned.a3m"
        sto_path="${OUTDIR}/${db}.sto"

        $PERL_PATH "${SCRIPT_DIR}/reformat.pl" a3m sto "$a3m_path" "$sto_path"
        python3 "${SCRIPT_DIR}/fix_reformat_sto.py" "$sto_path" "$sto_path"
    done

    echo "    [✓] STOs generated for $CLUSTER_NAME"
done

echo "[✓] All clusters processed and AF2 input-ready"
