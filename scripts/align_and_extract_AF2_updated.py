#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import numpy as np
import mdtraj as md
from typing import Tuple, List

def extract_ordered_ca_indices(traj: md.Trajectory, residue_range: Tuple[int, int]) -> List[int]:
    """Return atom indices of CA from specified inclusive residue range, ordered by resSeq."""
    start, end = residue_range
    ca = []
    for atom in traj.topology.atoms:
        if atom.name == "CA":
            rs = atom.residue.resSeq
            if start <= rs <= end:
                ca.append(atom.index)
    return ca

def main():
    p = argparse.ArgumentParser(description="Align PDBs to a reference on Cα and export flattened xyz to .npy")
    p.add_argument("ref_pdb", type=str, help="Reference PDB path used for alignment")
    p.add_argument("input_dir", type=str, help="Directory containing input .pdb files")
    p.add_argument("output_npy", type=str, help="Output .npy file")
    p.add_argument("--start", type=int, default=1, help="Start residue number (inclusive)")
    p.add_argument("--end", type=int, default=162, help="End residue number (inclusive)")
    p.add_argument("--pattern", type=str, default="*.pdb", help="Glob to select input files")
    p.add_argument("--labels_csv", type=str, default=None, help="Optional CSV to write row→filename mapping")
    p.add_argument("--strict", action="store_true", help="Abort if any structure is incomplete")
    args = p.parse_args()

    residue_range = (args.start, args.end)
    expected_len = args.end - args.start + 1

    print(f"[+] Reference: {args.ref_pdb}")
    ref = md.load(args.ref_pdb)
    ref_idx = extract_ordered_ca_indices(ref, residue_range)
    if len(ref_idx) != expected_len:
        raise ValueError(f"Reference missing CA atoms in {residue_range}. Found {len(ref_idx)} expected {expected_len}.")
    ref_ca = ref.atom_slice(ref_idx)

    files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    print(f"[+] Found {len(files)} PDB files in {args.input_dir} matching {args.pattern}")
    if not files:
        raise SystemExit("No input PDBs found.")

    rows = []
    kept = []
    for fp in files:
        try:
            traj = md.load(fp)
        except Exception as e:
            msg = f"[-] Skipping {os.path.basename(fp)} — failed to load: {e}"
            if args.strict:
                raise RuntimeError(msg)
            print(msg)
            continue

        idx = extract_ordered_ca_indices(traj, residue_range)
        if len(idx) != expected_len:
            msg = f"[-] Skipping {os.path.basename(fp)} — incomplete backbone ({len(idx)}/{expected_len} Cα)."
            if args.strict:
                raise RuntimeError(msg)
            print(msg)
            continue

        ca = traj.atom_slice(idx)
        ca.superpose(ref_ca)
        rows.append(ca.xyz[0].flatten())
        kept.append(os.path.basename(fp))

    if not rows:
        raise SystemExit("No usable structures. Nothing written.")

    arr = np.vstack(rows)
    np.save(args.output_npy, arr)
    print(f"[✓] Saved {arr.shape[0]} aligned structures to {args.output_npy} with shape {arr.shape}")

    if args.labels_csv:
        import csv
        with open(args.labels_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["row", "filename", "cluster"])
            for i, name in enumerate(kept):
                # best effort cluster id extraction, e.g. T4-Lysozyme_U100-004_ranked_3.pdb → U100-004
                cluster = ""
                base = os.path.splitext(name)[0]
                parts = [p for p in base.split("_") if p.startswith("U") and "-" in p]
                if parts:
                    cluster = parts[0]
                w.writerow([i, name, cluster])
        print(f"[✓] Wrote labels to {args.labels_csv}")

if __name__ == "__main__":
    # Maintain backwards compatible CLI
    if len(sys.argv) == 4 and sys.argv[1].endswith(".pdb"):
        # Original usage: python align_and_extract_backbones.py <ref_pdb> <input_dir> <output_npy>
        sys.argv = [sys.argv[0],
                    sys.argv[1], sys.argv[2], sys.argv[3],
                    "--start", "1", "--end", "162"]
    main()
