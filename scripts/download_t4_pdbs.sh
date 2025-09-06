#!/bin/bash
# === download_t4_pdbs.sh ===
# Downloads PDB structures from list into T4-Lysozyme directory

ID_LIST="/data/home/ap24012/structure_refs/T4-Lysozyme/t4_pdb_ids.txt"
OUT_DIR="/data/home/ap24012/structure_refs/T4-Lysozyme"

mkdir -p "$OUT_DIR"

while read -r pdbid; do
    pdbid_lower=$(echo "$pdbid" | tr '[:upper:]' '[:lower:]')
    url="https://files.rcsb.org/download/${pdbid_lower}.pdb"
    out_file="${OUT_DIR}/${pdbid}.pdb"
    
    if [[ ! -f "$out_file" ]]; then
        echo "[*] Downloading $pdbid from $url"
        wget -q -O "$out_file" "$url"
        if [[ $? -ne 0 ]]; then
            echo "[!] Failed to download $pdbid"
            rm -f "$out_file"
        fi
    else
        echo "[=] Skipping $pdbid (already exists)"
    fi
done < "$ID_LIST"

echo "[✓] Download complete. Files saved in: $OUT_DIR"
