#!/usr/bin/env python3

import os
import sys
import numpy as np
import mdtraj as md
import glob

def extract_ordered_ca_indices(traj, residue_range):
    """Return atom indices of CAs from specified residue range, in order."""
    ca_indices = []
    for atom in traj.topology.atoms:
        if atom.name == "CA" and residue_range[0] <= atom.residue.resSeq <= residue_range[1]:
            ca_indices.append(atom.index)
    return ca_indices

def align_and_extract(folder, ref_path, residue_range=(1,162), output_npy="backbone_aligned.npy"):
    print(f"[+] Reference: {ref_path}")
    print(f"[+] Input folder: {folder}")
    print(f"[+] Output file: {output_npy}")

    # Load and extract reference
    ref_traj = md.load(ref_path)
    ref_indices = extract_ordered_ca_indices(ref_traj, residue_range)
    ref_ca = ref_traj.atom_slice(ref_indices)

    files = sorted(glob.glob(os.path.join(folder, "*.pdb")))
    print(f"[+] Found {len(files)} PDB files.")

    data = []
    used = 0
    for pdb in files:
        traj = md.load(pdb)
        indices = extract_ordered_ca_indices(traj, residue_range)

        # Skip if missing atoms
        if len(indices) != (residue_range[1] - residue_range[0] + 1):
            print(f"[-] Skipping {os.path.basename(pdb)} — incomplete backbone.")
            continue

        ca = traj.atom_slice(indices)
        ca.superpose(ref_ca)
        data.append(ca.xyz[0].flatten())
        used += 1

    arr = np.array(data)
    np.save(output_npy, arr)
    print(f"[✓] Saved {used} aligned structures to {output_npy} with shape {arr.shape}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python align_and_extract_backbones.py <ref_pdb> <input_dir> <output_npy>")
        sys.exit(1)

    ref_path = sys.argv[1]
    input_dir = sys.argv[2]
    output_file = sys.argv[3]

    align_and_extract(input_dir, ref_path, (1,162), output_file)
