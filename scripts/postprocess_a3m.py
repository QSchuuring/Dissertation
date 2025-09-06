#!/usr/bin/env python3
import sys
from collections import OrderedDict

def read_a3m(path):
    sequences = OrderedDict()
    with open(path) as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_id = line
                if current_id in sequences:
                    # Duplicate ID, skip this one
                    current_id = None
                else:
                    sequences[current_id] = ""
            elif current_id:
                sequences[current_id] += line
    return sequences

def pad_sequences(sequences):
    max_len = max(len(seq) for seq in sequences.values())
    return {k: v.ljust(max_len, '-') for k, v in sequences.items()}

def write_a3m(sequences, path):
    with open(path, "w") as f:
        for header, seq in sequences.items():
            f.write(f"{header}\n{seq}\n")

if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    seqs = read_a3m(input_path)
    seqs = pad_sequences(seqs)
    write_a3m(seqs, output_path)
