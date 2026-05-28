#!/usr/bin/env python3
import sys,os
import re
import textwrap # for wrapping sequence lines
from datetime import datetime
from collections import defaultdict
from urllib.parse import unquote # Geneious output is URL-quoted
from typing import TextIO,Generator
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

# richer messages are in vogue
GREEN = '\033[92m'
RED = '\033[91m'
GREY = '\033[90m'
BOLD = '\033[1m'
RESET = '\033[0m'



metadata = defaultdict(dict) # map by seq_id: { header_field: value }
def argmax(l):
    max_el = max(l)
    return l.index(max_el)

pat = re.compile(r'N+')

def find_inner_chunk_no_N(seq):
    # split on Ns, take largest as trimmed sequence
    chunks = re.split(pat, seq)
    largest_chunk_ix = argmax([len(chunk) for chunk in chunks])
    largest_chunk = chunks[ largest_chunk_ix ]

    # validate
    flanks = seq.split(largest_chunk)
    # this should yield the flanks of leftFlank, largest_chunk, rightFlank
    # but if the largest chunk isn't unique, the result will be larger than 
    # 2 and this might not work
    if len(flanks) != 2:
        raise ValueError

    # return bounds
    s = len(flanks[0])
    e = s + len(largest_chunk)
    # core sequence will be seq[s:e]
    return s, e

def examine_sequence(seq:str, start:int, end:int, phase:int) -> tuple:
    if seq.find('N') >= 0:
        trimmed_range = find_inner_chunk_no_N(seq)
        # adjust coordinates
        # start = start - trimmed_range[0]
        # end = end - trimmed_range[0]
        # seq = seq[trimmed_range[0]:trimmed_range[1]]

    if end > len(seq):
        raise ValueError
    
    cds = seq[(start-1)+phase:end]  # convert to 0-based
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    has_start = cds[:3].upper() == 'ATG'
    has_stop  = cds[-3:].upper() in stop_codons
    is_5prime_partial = not has_start
    is_3prime_partial = not has_stop
    
    # codon_start: how many bases into the first codon are we?
    # If 5' partial and length % 3 != 0, you need to figure out the offset
    # Usually set to 1 unless you have evidence otherwise
    codon_start = phase + 1
    #if is_5prime_partial and not is_3prime_partial: # here is a situation where the STOP codon exists but with no START codon,
                                                    # so the 1-based frameshift is gotten from the modulus
        #codon_start = (len(seq) % 3) + 1
    
    return is_5prime_partial, is_3prime_partial, codon_start

def write_tbl_entry(genbank_id:str, partial3prime: bool, partial5prime: bool, start: int, end: int, codon_start: int, gene_name:str, tbl_out: TextIO) -> None:
    s = f"<{start}" if partial5prime else str(start)
    e = f">{end}"   if partial3prime else str(end)
    print(f">Feature {genbank_id}", file=tbl_out)
    print(s, e, "CDS", sep="\t", file=tbl_out)
    print("\t" * 3 + "product", gene_name, sep="\t", file=tbl_out)
    print("\t" * 3 + "codon_start", codon_start, sep="\t", file=tbl_out)
    print("\t" * 3 + "transl_table", 1, sep="\t", file=tbl_out)

import glob

def main():

    ### hard-code paths- this is not a generalizable script
    basedir= os.path.dirname(sys.argv[0]) # where are we?

    # get prefix (like "current/A2")
    if len(sys.argv) < 2:
        user_prefix = "current/send_to_david/FABADRU_SA"
    else:
        user_prefix = sys.argv[1]

    out_dir = f"{user_prefix}_submission_files"
    if os.path.exists(out_dir):
        print(f"Output director {out_dir} exists")
    else:
        os.mkdir(out_dir)
        print(f"Made output directory {out_dir}")

    gff_file = glob.glob(f"{user_prefix}*.gff")[0]
    tsv_file = glob.glob(f"{user_prefix}*.txt")[0]
    fasta_file = glob.glob(f"{user_prefix}*.fasta")[0]

    outfile_prefix = os.path.basename(user_prefix)
    tbl_outpath = f"{out_dir}/{outfile_prefix}.tbl"
    fsa_outpath = f"{out_dir}/{outfile_prefix}.fsa"
    txt_outpath = f"{out_dir}/{outfile_prefix}.txt"
    gff_outpath = f"{out_dir}/{outfile_prefix}.gff"

    GFF = tillie.read_GFF(gff_file)

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


    with open(fsa_outpath, "w") as fsa_out, open(tbl_outpath, "w") as tbl_out:

        for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
            seq_id = record.description.strip()
            parsed_seq_id = tillie.parse_seq_header(seq_id)
            genbank_id = tillie.canonical_seq_name(seq_id)
            matched_row_in_metadata = metadata_df[metadata_df['id_for_genbank'] == genbank_id]

            # crash if no match found
            if len(matched_row_in_metadata) != 1:
                if len(matched_row_in_metadata) == 0:
                    print(f"Error, sequence {seq_id} matched no rows in metadata looking for strain = {parsed_seq_id['strain']} AND prot = {parsed_seq_id['prot']}", file=sys.stderr)   
                # else:
                    # match prot to gene_name to get the right GFF row (there are erroneous annotations)
                raise ValueError

            # determine CDS features
            gff = GFF[ GFF['seqid'] == genbank_id]
            start = gff['start'].item()
            end = gff['end'].item()
            gene_name = gff['gene_name'].item()
            phase = 0
            try:
                phase = int(gff['phase'].item())
            except:
                pass
            
            #codon_start = int(gff['phase'].item()) + 1
            partial5prime, partial3prime, codon_start = examine_sequence(str(record.seq), start, end, phase)
            # If we found a start codon, override phase to codon_start=1
            if not partial5prime:
                codon_start = 1

            # create a new entry for the .tbl file
            write_tbl_entry(genbank_id, partial3prime, partial5prime, start, end, codon_start, gene_name, tbl_out)

            # the sequence header has more specific host values, use if available
            if 'host' in parsed_seq_id:
                matched_row_in_metadata['host'] = parsed_seq_id['host']

            # save new row and remove from original
            processed_df = pd.concat([processed_df, matched_row_in_metadata], axis=0)
            metadata_df = metadata_df.drop(matched_row_in_metadata.index)

            matched_seq_id = matched_row_in_metadata['seq_ID'].item()
            annotations = []
            export_fields_keys = ["genotype", "host", "strain", "collected_by", "geo_loc_name", "Collection_date", "note", "isolation_source"] # 'gene' in the sequence file is unsupported in future releases
        
            for k,v in matched_row_in_metadata.items():
                if k in export_fields_keys:
                    if k == 'genotype' and v.item() == 'NA': continue
                    if k == 'Collection_date' and v.item() == 'NA': continue

                else: 
                    continue

                val = str(v.item()).replace('"','')
                annotations.append( f"[{str(k).strip()}={val}]" )
            annotations.append( '[organism=Bluetongue virus]')
            # write annotated header and sequence
            print('>' + genbank_id, *annotations, file=fsa_out)
            print(str(record.seq), file=fsa_out)

    print(f"{GREEN}Processed {i + 1} sequences from {BOLD}{fasta_file}.{RESET}", file=sys.stderr)
    print("Done!", file=sys.stderr)
    print(f"Wrote {fsa_outpath}", file=sys.stderr)
    print(f"Wrote {tbl_outpath}", file=sys.stderr)

    # extra genbank steps for the metadata file
    processed_df.replace("NA", "missing: lab stock", inplace=True)
    processed_df.drop(['BTV', 'number_segment_sequences','product','seq_ID'], axis=1, inplace=True)
    processed_df.rename(columns={'id_for_genbank': 'seq_ID'}, inplace=True)
    processed_df['Organism'] = "Bluetongue virus"

    processed_df.to_csv(txt_outpath, sep="\t", index=False)


    with open(gff_outpath,"w") as gff_out:
        print("##gff-version 3", file=gff_out)
        GFF.to_csv(gff_out, sep="\t", index=False, mode="a", header= False)
    if len(metadata_df) > 0:
        print("left over in metadata:", file=sys.stderr)
        metadata_df.to_csv(sys.stderr, sep="\t", index=False)
    else:
        print("No unaccounted for metadata (good).", file=sys.stderr)

    print(f"Wrote {txt_outpath}", file=sys.stderr)

################


def iterate_and_replicate_rows(fields:dict[str,str]) -> Generator[dict, None, None]:
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