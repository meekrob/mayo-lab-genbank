#!/usr/bin/env python3
import sys
import pandas as pd
import re
import textwrap
from collections import defaultdict
from Bio import SeqIO
from datetime import datetime
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

metadata = defaultdict(lambda: {}) # map by seq_id: { header_field: value }

from genbank_vals import geoloc_countries, state_names_to_abbrev, host_lookup

def id_seqid_type(fields, btv_col):
    """
    fields - parsed from sequence header
    btv_col - the column index containing the "BTV" entry. -1 if missing.
    
    The input fasta headers seem to come in 4 different flavors depending on what data are available.
    Only "prot" and "strain" are common to all, but those two are sufficient to map to entries in the metadata file.

    """
    # most scant version of row. just FABRADU
    if btv_col == -1:
        prot, strain = fields
        return {
            'prot': prot,
            'strain': strain #.replace('_','') # needed in previous version of data
        }
    
    # VP3     Clinical_16     BTV6    mule_deer       5.5yo   female  2021.10.19      2/7
    if btv_col == 2 and len(fields) == 7:
        prot, strain, btv, host, age, sex, date = fields
        return { 'prot': prot,
            'strain': strain,
            'btv': btv,
            'host': host,
            'age': age,
            'sex': sex,
            'date': date.replace('.', '-')}

    # VP3     BTV17   AA17AAAFAAAC    CA77_Caramel    Bovine  Kern_California 2012-10-01      1/7    
    if btv_col == 1 and len(fields) == 7:
        
        return {
            'prot': fields[0],
            'btv': fields[1],
            'genotype': fields[2],
            'strain': fields[3],
            'host': fields[4],
            'location': fields[5],
            'date': fields[6]
        }
    
    # VP3     BTV22   AG22AAAAAAAA    FABADRU000164   Trinidad        1/5
    if btv_col == 1 and len(fields) == 5:
        prot, btv, genotype, strain, location = fields
        return {
            'prot': prot,
            'btv': btv,
            'genotype': genotype,
            'strain': strain,
            'location': location
        }
    
    return { 'error': ';'.join(fields)}

def parse_seq_header(header_line):
    fields = re.split('[/|]', header_line.strip())
    found_BTV = False
    for i, f in enumerate(fields):
        if f.startswith('BTV'):
            found_BTV = True
            parsed = id_seqid_type(fields, i)
            return parsed
            
    if not found_BTV:
        parsed = id_seqid_type(fields, -1)
        return parsed

    print("Didn't recognize format of seq header:", header_line, file=sys.stderr)    
    raise ValueError

def print_parsed_header(parsed_dict):
    for k,v in parsed_dict.items():
        print(f"{k}={v}", end=' ')

    print()


def convert_date_and_location(row):

    # combine columns to make the geo location
    if row['country'] == 'United States':
        row['country'] = 'USA'

    elif row['country'] == 'Trinidad':
        row['country'] = 'Trinidad and Tobago'

    elif row['country'] not in geoloc_countries:
        print(f"Country label {row['country']} not in list.", file=sys.stderr)
        sys.exit(1)

    row['state'] = state_names_to_abbrev[ row['state'] ]

    row['geo_loc_name'] = row['country']
    if row['state']:
        if row['county'] == 'NA':
            county = ''
        else:
            county = row['county'] + ' County, '
        row['geo_loc_name'] += ": " + county + row['state']

    del row['state']
    del row['county']
    del row['country']
    del row['continent']

    # combine columns to make the date
    if "NA" in [row['year'],row['month'],row['day']]:
        date = "NA"
    else:
        date = datetime(int(row['year']), int(row['month']), int(row['day'])).strftime("%d-%b-%Y")
    row['Collection_date'] = date
    del row['year']
    del row['month']
    del row['day']

    return row



platform_mapping = {
    "Nanopore Sequencing": "nanopore",
    "Illumina Sequencing": "illumina"
}

isolation_source_map = {
    "cell culture isolate": "cell culture",
    "Lung": "lung",
    "Blood": "blood",
    "Spleen": "spleen",
    "Lymph node": "lymph node",
    "whole blood": "whole blood",
    '"Pool: Lung, Spleen, & Lymph"': "pooled_sample",
    "Pool: Lung, Spleen, & Lymph": "pooled sample",
    "Pool: Lung & Spleen": "pooled sample"
}

def validate_genbank_fields(fields):
    if 'sequencing_platform' in fields:
        fields['note'] = f"sequenced with {platform_mapping[ fields['sequencing_platform'] ]}"
        del fields['sequencing_platform']

    if 'source' in fields:
        fields['isolation_source'] = isolation_source_map[ fields['source'] ]
        if fields['source'].startswith("Pool"):
            source_note = fields['source'].lower().replace(' &', ',')
            if 'note' in fields:
                fields['note'] += '; ' + "pooled from lung, spleen, lymph"
            else:
                fields['note'] = "pooled from lung, spleen, lymph"

        del fields['source']

    return fields

def iterate_and_replicate_rows(fields):
    """
    iterate_and_replicate_rows - take a metadata row and repeat it for each of the gene products specified by the "number_segment_sequences" column
    """
    fields = convert_date_and_location(fields)
    fields = validate_genbank_fields(fields)

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

                # yield from here
                #print(f"Expanding {metadata[seq_id]['seq_ID']}:", file=sys.stderr)
                for processed_row in iterate_and_replicate_rows(metadata[seq_id]):
                    metadata[processed_row['seq_ID']] = processed_row
                    #print(f"{records}\t{processed_row['seq_ID']}:", ";".join(metadata[processed_row['seq_ID']].values()))
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
processed_df = pd.DataFrame() # transfer rows here when processed
track_seq_id_unique = defaultdict(lambda: [])

for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
    seq_id = record.description.strip()
    parsed_seq_id = parse_seq_header(seq_id)
    matched_row_in_metadata = metadata_df[(metadata_df['product'] == parsed_seq_id['prot']) & (metadata_df['strain'] == parsed_seq_id['strain'])]
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

    #annotations.append(f"[orig_seq_name = {seq_id}]")
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