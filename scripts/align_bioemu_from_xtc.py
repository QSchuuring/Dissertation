import mdtraj as md
import numpy as np
import sys

if len(sys.argv) != 5:
    print("Usage: python align_bioemu_from_xtc.py <ref_pdb> <topology.pdb> <trajectory.xtc> <output.npy>")
    sys.exit(1)

ref_path = sys.argv[1]
topo_path = sys.argv[2]
traj_path = sys.argv[3]
out_path = sys.argv[4]

# Residues 1–162 (inclusive), which means resid 0–161 in MDTraj
res_start, res_end = 0, 161
expected_n_residues = res_end - res_start + 1  # 162 residues

# Load reference structure and extract CAs
ref = md.load(ref_path)
ref_ca = ref.atom_slice(ref.topology.select(f"name CA and resid >= {res_start} and resid <= {res_end}"))

# Load trajectory
traj = md.load_xtc(traj_path, top=topo_path)

# Extract CA atoms
ca_indices = traj.topology.select(f"name CA and resid >= {res_start} and resid <= {res_end}")
print(f"[INFO] Found {len(ca_indices)} CA indices")

traj_ca = traj.atom_slice(ca_indices)

# Check that each frame has exactly 162 atoms
if traj_ca.n_atoms != expected_n_residues:
    print(f"[-] Mismatch: found {traj_ca.n_atoms} CAs, expected {expected_n_residues}.")
    sys.exit(1)

# Align each frame to reference
traj_ca.superpose(ref_ca)

# Flatten and save
coords = traj_ca.xyz.reshape(traj_ca.n_frames, -1)
np.save(out_path, coords)
print(f"[✓] Saved BioEmu aligned trajectory → {out_path} with shape {coords.shape}")
