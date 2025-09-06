#!/bin/bash
# submit_af2_pipeline.sh
# === Dynamic batch submission of AlphaFold2 jobs for clustered MSAs ===

# === CONFIGURATION ===
PROT="$1"  # e.g. lysozyme
FASTA=/data/home/ap24012/alphafold_input/${PROT}/fastas/${PROT}.fasta
CLUSTER_DIR=/data/home/ap24012/scratch_preserved/af2_runs/${PROT}/clusters
MSA_SOURCE_DIR=/data/home/ap24012/scratch_preserved/af2_runs/${PROT}  # where each cluster folder will be created
RUN_SCRIPT=/data/home/ap24012/alphafold_scripts/run.sh
TEMP_BASE=/data/home/ap24012/scratch_preserved/af2_output
FINAL_BASE=/data/home/ap24012/scratch_preserved/af2_results

# === PROCESS EACH CLUSTERED MSA ===
for size in $(seq 12 2 30); do
    for A3M in ${CLUSTER_DIR}/${PROT}_U${size}-*.a3m; do
        [ -e "$A3M" ] || continue   # skip if none exist
        CLUSTER=$(basename "$A3M" .a3m)

        # Directories
        CLUSTER_DIR_OUT="${MSA_SOURCE_DIR}/${CLUSTER}"
        TEMP_DIR="${TEMP_BASE}/${CLUSTER}"
        FINAL_DIR="${FINAL_BASE}/${CLUSTER}"

        ALIGNED="${CLUSTER_DIR_OUT}/${CLUSTER}_aligned.a3m"
        STO1="${CLUSTER_DIR_OUT}/uniref90_hits.sto"
        STO2="${CLUSTER_DIR_OUT}/small_bfd_hits.sto"
        STO3="${CLUSTER_DIR_OUT}/mgnify_hits.sto"

        echo "[*] Submitting AF2 job for ${CLUSTER}"

        qsub -N AF2_${CLUSTER} <<EOF
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=1:00:00
#$ -l h_vmem=11G
#$ -l gpu=1

# === Environment Setup ===
export OUTPUT=${TEMP_DIR}
module load alphafold2

mkdir -p \${OUTPUT}/${PROT}/msas
cp "${STO1}" \${OUTPUT}/${PROT}/msas/uniref90_hits.sto
cp "${STO2}" \${OUTPUT}/${PROT}/msas/small_bfd_hits.sto
cp "${STO3}" \${OUTPUT}/${PROT}/msas/mgnify_hits.sto

alphafold ${RUN_SCRIPT} \\
  --fasta_paths=${FASTA} \\
  --use_precomputed_msas=True \\
  > \${OUTPUT}/${CLUSTER}_alphafold.log 2>&1

mkdir -p ${FINAL_DIR}
cp -r \${OUTPUT}/* ${FINAL_DIR}/
rm -rf \${OUTPUT}
EOF

    done
done

echo "[✓] All clustered AF2 jobs submitted for ${PROT}."
