#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) University of Strathclyde 2021
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
"""seqid_to_species_name.py

This script is an ad hoc script to make production of phylogenetic
trees with interpretable labels a little easier.

Ordinarily, the first element in a FASTA header is a unique sequence
ID. This is done for good reasons, but as it is used to label trees
produced from a sequence file, can result in a bit of confusion when
trying to interpret those trees.

This script will take a FASTA file and a corresponding .csv file
formatted as follows:

    seqid,species_name
    [seqID1],[species_name1]
    [seqID2],[species_name2]
    ...

and use the .csv file to replace the seqID with the species_name in a
new output FASTA file.

The old sequence ID, and the description string, are preserved as the
description string in the new FASTA file.

As sequence IDs should not contain whitespace (whitespace indicates
the start of the description string in a FASTA header), any whitespace
in the proposed species_name is replaced with an underscore (_).
"""

import csv
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

from Bio import SeqIO


def parse_args():
    parser = ArgumentParser(
        prog="seqid_to_species_name.py", formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--fasta",
        dest="fasta_path",
        action="store",
        default=None,
        type=Path,
        help="Path to input FASTA file",
    )
    parser.add_argument(
        "--csv",
        dest="csv_path",
        action="store",
        default=None,
        type=Path,
        help="Path to input .csv file",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile_path",
        action="store",
        default=Path("out.fasta"),
        type=Path,
        help="Path to output FASTA file",
    )
    return parser.parse_args(sys.argv[1:])


def parse_csv(fpath):
    print(f"Reading new names from {fpath}")

    # Holds the new names as a dictionary, keyed by old name
    new_names = {}

    # Parse the CSV file
    with fpath.open("r") as ifh:
        reader = csv.reader(ifh)
        next(reader)  # skip the header
        for row in reader:
            new_names[row[0]] = row[1].replace(" ", "_")

    return new_names


def update_sequences(records, new_names):

    for record in records:

        # Ailsa has split the sequence ID before the sequence range
        # so we need to remove the sequence range
        seqid = record.id.split("/")[0]

        if seqid in new_names:
            record.description = f"{record.id} {record.description}"
            record.id = new_names[seqid]
        yield record


def main():
    # Get paths to input/output files
    args = parse_args()

    # Parse the CSV file to obtain new sequence names
    # as a dictionary, keyed by old name
    new_names = parse_csv(args.csv_path)

    # Parse the FASTA file as a generator
    print(f"Reading sequences from {args.fasta_path}")
    records = SeqIO.parse(args.fasta_path, "fasta")

    # Process the FASTA records, replacing the sequence ID
    # with the new name (receives a generator)
    print(f"Updating sequence IDs")
    new_records = update_sequences(records, new_names)

    # Write the new FASTA file
    print(f"Writing new FASTA file to {args.outfile_path}")
    SeqIO.write(new_records, args.outfile_path, "fasta")


if __name__ == "__main__":
    main()
