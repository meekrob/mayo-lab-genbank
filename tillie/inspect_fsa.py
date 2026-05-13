#!/usr/bin/env python3
"""
inspect_fsa.py - get sequences out of the formatted fsa file to help debug
"""
import sys
import tillie_data as tillie
from Bio import SeqIO

fasta_file = "/Users/david/work/mayo_lab_sequence_submission/tillie/additional_oc_sequences_for_submission.fsa"

fsa = {}
for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
    seq_id:str = record.id.strip()
    fsa[seq_id] = record

pass