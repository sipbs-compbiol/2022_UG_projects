#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path

import sys

# Catch filename from command-line
fname = Path(sys.argv[1])

# Open input file
with open(fname) as handle:
    records = list(SeqIO.parse(handle, "fasta"))

# Fix headers
for record in records:
    geneid = record.description.split("GN=")[-1].split(" ")[0]
    record.id = geneid

# Rename output file
new_stem = fname.stem + "_fixed"
with open(new_stem + ".fasta", "w") as handle:
    SeqIO.write(records, handle, "fasta")