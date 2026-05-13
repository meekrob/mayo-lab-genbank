#!/usr/bin/env python3
"""
summarize_fasta.py -
Get information about the sequences, whether the canonical naming is unique, etc.
"""
import sys
import tillie_data as tillie
from Bio import SeqIO
from collections import Counter

gff_file = "/Users/david/work/mayo_lab_sequence_submission/tillie/new_files/TD_n1257_GenBank.gff"
fasta_file = "/Users/david/work/mayo_lab_sequence_submission/tillie/new_files/TD_n1257_GenBank.fasta"

GFF = tillie.read_GFF(gff_file)

seqs_by_header = Counter()
seqs_by_canonical = Counter()

all_seqs = {}

for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
    seq_id:str = record.description.strip()
    parsed_seq_id:dict = tillie.parse_seq_header(seq_id)
    genbank_id:str = tillie.canonical_seq_name(seq_id)

    seqs_by_header[seq_id] += 1
    seqs_by_canonical[genbank_id] += 1
    all_seqs[genbank_id] = record


duplicate_headers = [ (k,count) for k,count in seqs_by_header.items() if count > 1]
duplicate_canonical = [ (k,count) for k,count in seqs_by_canonical.items() if count > 1]

print(len(duplicate_headers))
print(len(duplicate_canonical))

