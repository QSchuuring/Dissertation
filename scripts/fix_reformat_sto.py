fix_reformat_sto.py#!/usr/bin/env python

"""
Move #=GC line (if present) to just before the final // line in a Stockholm (.sto) file.

Usage:
python fix_reformat_sto.py <input.sto> <output.sto>

This script ensures that AlphaFold-compatible Stockholm files have any #=GC lines
placed just before the terminal // line, as required by the format.
"""

import os, sys

def main():
    argv = sys.argv
    input_filename = argv[1]
    output_filename = argv[2]

    # Read input .sto file
    with open(input_filename, 'r') as infile:
        infile_lines = infile.readlines()

    GCline = None  # Initialize in case no #=GC line exists

    # Write output, moving #=GC line (if any) to before //
    with open(output_filename, 'w') as output_file:
        for line in infile_lines:
            if line.startswith("#=GC"):
                GCline = line  # Save the GC annotation line for later
            elif line.startswith("//"):
                if GCline:
                    output_file.write(GCline)  # Insert saved GC line before //
                output_file.write(line)  # Write end-of-file line
            else:
                output_file.write(line)  # Write all other lines normally

if __name__ == "__main__":
    sys.exit(main())
