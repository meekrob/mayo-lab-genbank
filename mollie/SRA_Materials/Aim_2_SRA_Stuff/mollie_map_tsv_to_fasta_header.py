#!/usr/bin/env python3
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import mollie_data as mollie

USAGE = f"""
 This script is designed to help create FASTA files in the right format for Genbank submission.

 It creates a fasta "defline" in the appropriate format using information in a tab-delimited
 sequence modifier table.

 It expects a FASTA file to already exist with sequence names matching those in the first column of the 
 sequence modifier table.

 It expects the sequence modifier table to be tab-delimited and for the first row to contain column
 headers that are legitimate Genbank sequence modifiers.

 For info about sequence modifier tables, see:
      http://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html

 For info about FASTA deflines, see: 
      http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/#fsa

 For info about the genbank batch submission process using the tbl2asn utility, see:
      http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/#


  Usage: {sys.argv[0]} [-h] <sequence_modifier_table_file> <fasta_file>

         -h   print this help message

         Reads fasta file from stdin and writes modified fasta to stdout
 
   
  Mark Stenglein 4/21/2016
  Adapted to python by David King 5/4/2026 (original script didn't work on Molly Burton's files for some reason).
"""

if len(sys.argv) < 2 or sys.argv[1] == '-h':
    print(USAGE)
    sys.exit(0)


# richer messages are in vogue
GREEN = '\033[92m'
RED = '\033[91m'
GREY = '\033[90m'
BOLD = '\033[1m'
RESET = '\033[0m'

#tsv_file = sys.argv[1]
import os
tsv_file = os.path.dirname(sys.argv[0])
#fasta_file = sys.argv[2]

metadata = defaultdict(lambda: {}) # map by seq_id: { header_field: value }

records = 0
with open(tsv_file) as tsv_fh:
    header = []
    for line_i, line in enumerate(tsv_fh):
        if line.startswith('#'): continue
    
        fields = line.strip().split("\t")

        # first line encountered
        if not header:
            header = fields[1:]
            print(f"{GREEN}{tsv_file}: Encountered header at index: {line_i}{RESET}", file= sys.stderr)
            continue

        # all else a data line
        if len(line.strip()) > 0:
            try:
                seq_id = fields[0]
                for k,v in zip(header, fields[1:]):
                    metadata[seq_id][k] = v

                records += 1

            except IndexError:
                print(f"Error parsing data line from {tsv_file}.", file=sys.stderr)
                print(f"line: >{line}<", file=sys.stderr)
                print("Fields:",  fields, file=sys.stderr)
                raise
        else:
            print(f"{GREY}{tsv_file}: Encountered blank line at index: {line_i}{RESET}", file= sys.stderr)

    print(f"{GREEN}Processed {records} records from {line_i+1} lines from {BOLD}{tsv_file}.{RESET}",  file=sys.stderr)
 

for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
    seq_id = record.description.strip()
    if seq_id in metadata:
        print('>' + seq_id.replace(' ', '_'), end=' ')
        for h in header:
            print(f"[{h}={metadata[seq_id][h]}]", end=" ")
        print()
    else:
        print(f"Warning: {seq_id} not in data file", file=sys.stderr)
        continue

    print(str(record.seq))

print(f"{GREEN}Processed {i + 1} sequences from {BOLD}{fasta_file}.{RESET}", file=sys.stderr)
print("Done!", file=sys.stderr)
