#!/usr/bin/env python3
import sys
import re
import textwrap # for wrapping sequence lines
from datetime import datetime
from collections import defaultdict
from urllib.parse import unquote # Geneious output is URL-quoted
import pandas as pd
from Bio import SeqIO

import tillie_data as tillie # to parse and canonicalize naming conventions she used

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
  Adapted further by David King 5/7/2026 to replicate rows, adding protein names to seq ID
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

tsv_file = sys.argv[1]
fasta_file = sys.argv[2]
gff_file = sys.argv[3]

metadata = defaultdict(dict) # map by seq_id: { header_field: value }


def main():

    # gff
    GFF = pd.read_csv(gff_file, sep="\t", header = None, comment='#')
    GFF.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    GFF['seqid'] = GFF['seqid'].apply(unquote)
    GFF['seqid'] = GFF['seqid'].apply(tillie.parse_seq_header)

    records = 0
    with open(tsv_file) as tsv_fh:
        header = []
        for line_i, line in enumerate(tsv_fh):
            if line.startswith('#'): continue
        
            fields = []
            for f in line.strip().split("\t"):
                if f == '': # blank cell means end of data in row
                    break
                else: 
                    fields.append(f.strip())

            # first line encountered
            if not header:
                header = fields
                print(f"{GREEN}{tsv_file}: Encountered header at index: {line_i}{RESET}", file= sys.stderr)
                continue

            # all else a data line
            if len(line.strip()) > 0:
                try:
                    seq_id = fields[0]
                    for k,v in zip(header, fields):
                        if v == '': break # for the rows with notes off to the right (cell is empty), end of data was reached
                        metadata[seq_id][k] = v

                    for processed_row in iterate_and_replicate_rows(metadata[seq_id]):
                        metadata[processed_row['seq_ID']] = processed_row
                        records += 1
                    del metadata[seq_id]

                except IndexError:
                    print(f"Error parsing data line from {tsv_file}.", file=sys.stderr)
                    print(f"line: >{line}<", file=sys.stderr)
                    print("Fields:",  fields, file=sys.stderr)
                    raise
            else:
                print(f"{GREY}{tsv_file}: Encountered blank line at index: {line_i}{RESET}", file= sys.stderr)

        print(f"{GREEN}Processed {records} records from {line_i+1} lines from {BOLD}{tsv_file}.{RESET}",  file=sys.stderr)
    
    metadata_df = pd.DataFrame(metadata.values())
    metadata_df['id_for_genbank'] = metadata_df['seq_ID'].apply(tillie.canonical_seq_name)
    processed_df = pd.DataFrame() # transfer rows here when processed
    track_seq_id_unique = defaultdict(list)

    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        seq_id = record.description.strip()
        parsed_seq_id = tillie.parse_seq_header(seq_id)
        genbank_id = tillie.canonical_seq_name(seq_id)
        #matched_row_in_metadata = metadata_df[(metadata_df['product'] == parsed_seq_id['prot']) & (metadata_df['strain'] == parsed_seq_id['strain'])]
        matched_row_in_metadata = metadata_df[metadata_df['id_for_genbank'] == genbank_id]
        if len(matched_row_in_metadata) != 1:
            if len(matched_row_in_metadata) == 0:
                print(f"Error, sequence {seq_id} matched no rows in metadata looking for strain = {parsed_seq_id['strain']} AND prot = {parsed_seq_id['prot']}", file=sys.stderr)   
            raise ValueError


        processed_df = pd.concat([processed_df, matched_row_in_metadata], axis=0)
        metadata_df = metadata_df.drop(matched_row_in_metadata.index)

        matched_seq_id = matched_row_in_metadata['seq_ID'].item()
        annotations = []
        for k,v in matched_row_in_metadata.items():
            if k == 'seq_ID': continue
            if k == 'prot': continue
            if k == 'genotype' and v.item() == 'NA': continue
            if k == 'product': 
                k = 'gene'
            if k == 'number_segment_sequences':
                continue
            if k == 'BTV':
                continue
            else:
                val = str(v.item()).replace('"','')
                annotations.append( f"[{str(k).strip()}={val}]" )

        # write annotated header and sequence
        product = matched_row_in_metadata['product'].item()
        btv = matched_row_in_metadata['BTV'].item()
        strain = matched_row_in_metadata['strain'].item()
        sanitized_seq_id = f"{product}_{btv}_{strain}"

        # the following is to make sure the sanitized_seq_id stays unique
        if sanitized_seq_id in track_seq_id_unique:
            print(f"error: {sanitized_seq_id} already used for {','.join(track_seq_id_unique[sanitized_seq_id])}", file = sys.stderr)
            sys.exit(1)
        
        track_seq_id_unique[sanitized_seq_id].append(matched_seq_id)

        # done, print new sequence info
        print('>' + sanitized_seq_id, *annotations)
        print(textwrap.fill(str(record.seq),width=150))

    print(f"{GREEN}Processed {i + 1} sequences from {BOLD}{fasta_file}.{RESET}", file=sys.stderr)
    print("Done!", file=sys.stderr)

    print("left over in metadata:", file=sys.stderr)
    metadata_df.to_csv(sys.stderr, sep="\t", index=False)

################


def iterate_and_replicate_rows(fields):
    """
    iterate_and_replicate_rows - take a metadata row and repeat it for each of the gene products specified by the "number_segment_sequences" column
    """
    fields = tillie.convert_date_and_location(fields)
    fields = tillie.validate_genbank_fields(fields)

    # need BTVXX out of the ID but it will have to be dropped before final printout
    fields['BTV'] = fields['seq_ID'].split('/')[0]

    if 'number_segment_sequences' in fields:
        
        try:
            number_segment_sequences = fields['number_segment_sequences']
            if int(number_segment_sequences) == 10:
                # replicate the row by prepending the following Protein abbreviation to the seq ID
                prots = [ 'VP1', 'VP2', 'VP3', 'VP4', 'NS1', 'VP5', 'VP7', 'NS2', 'VP6', 'NS3' ]
        except ValueError:
            num_seqs, seq_list = number_segment_sequences.strip('"').rstrip(')').split('(')
            prots = seq_list.split(',')
            if int(num_seqs) != len(prots):
                print(f"error, number of sequences named {num_seqs} does not match provided list {seq_list}", file=sys.stderr)
                sys.exit(1)
            
        for prot in prots: 
            f = fields.copy()
            f['seq_ID'] = prot.strip() + "/" + fields['seq_ID']
            f['product'] = prot.strip()
            yield f
            
    else:
        print(f"Error, field: 'number_of_sequences' not found in row", file=sys.stderr)
        raise IndexError


if __name__ == "__main__": main()