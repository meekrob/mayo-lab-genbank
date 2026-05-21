# tillie_data.py
import sys
import re
from datetime import datetime
from genbank_vals import geoloc_countries, state_names_to_abbrev, host_lookup
import pandas as pd
from Bio import SeqIO
from urllib.parse import unquote # Geneious output is URL-quoted

def get_gene_name_from_spreadsheet_sampleID(sampleID:str) -> str:
    if sampleID.startswith('Colorado deer'): return 'RdRP'
    if sampleID.startswith('Bovine hepacivirus'): return "Polyprotein"
    if sampleID.startswith('Totivirus nijyuroku'): return 'pol'
    if sampleID.endswith("polymerase"): return "RNA-dependent"
    if sampleID.endswith("hypothetical protein"): return "hypothetical protein"
    if sampleID.endswith("capsid"): return "putative capsid"

    return sampleID.strip().split()[-1]

def get_gene_name(attributes:str) -> str | None:
    match = re.search(r'Name=([^;]+)', attributes) # assume gene_name is indentified by Name=... with an optional ';' 
    if match and match.group(1).strip():
        gene_name = match.group(1).replace(' CDS', '')
        return gene_name
    return None

def get_seq_id_prot(orig_seq_id:str) -> str:
    # this is a special function to get just the gene name as specified by the seq ID in the GFF file. 
    # It needs to be called separately in order to contrast it with the gene_name parsed from the attribute column
    parsed = parse_seq_header(orig_seq_id, "")
    return parsed['prot']

def read_GFF(gff_file) -> pd.DataFrame:
    """
    read_GFF -
    Use pandas dataframes to filter and format the Geneious-output annotations
    1. Only CDS annotations are used
    2. Unique records must be enforced so that the mapping to the fasta file is 1-to-1
    3. The function canonical_seq_name() detects the original seqname form and make them consistent across input files
    """
    
    GFF = pd.read_csv(gff_file, sep="\t", header = None, comment='#')
    GFF.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    # filter out Geneious "Editing History..." feature rows
    # and empty Name= attribute rows (present in Tillie's data, idk about Mollie's).
    GFF = GFF[ (GFF['type'] == 'CDS') & 
              (~GFF['attributes'].str.contains(r'Name=$', regex=True))].copy() 

    # extract gene name from attributes, but fail if not found
    GFF['gene_name'] = GFF['attributes'].apply(get_gene_name)
    missing_gene_name = GFF[GFF['gene_name'].isna()]
    if not missing_gene_name.empty:
        print("CDS rows with unparseable gene names:", file=sys.stderr)
        print(missing_gene_name[['seqid', 'attributes']], file=sys.stderr)
        raise ValueError(f"{len(missing_gene_name)} CDS rows failed gene name extraction")
    
    # parse seqid into our canonical one
    GFF['seqid'] = GFF['seqid'].apply(unquote)
    GFF['orig_id'] = GFF['seqid']
    GFF['seqid'] = GFF.apply(
        lambda row: canonical_seq_name(row['seqid'], row['gene_name']),
        axis = 1
    )

    return GFF

def canonical_seq_name_from_parsed(parsed:dict) -> str:
    canonical_name = parsed['canonical']
    if len(canonical_name) > 40:
        print(f"name is too long: {canonical_name}: {len(canonical_name)}", file=sys.stderr)
        raise ValueError
    return canonical_name

def canonical_seq_name(header_line:str, gene_name:str) -> str:
    parsed = parse_seq_header(header_line, gene_name)
    if gene_name is not None:
        parsed['prot'] = gene_name
    return canonical_seq_name_from_parsed(parsed)

def canonical_host(host_field) -> str:
    host = host_field.lower()
    if host in host_lookup:
        return host_lookup[host]
    
    return host

def get_virus_prefix(term):
    # these have the abbreviations parsed
    if term.startswith('Skunk'): return ''
    if term.startswith('Black'): return ''
    if term.startswith('Bluetongue'): return ''
    # these lose the abbreviation during parsing
    if term.startswith('Horsetooth deer virus'): return "HDV_"
    if term.startswith('Totivirus nijyuroku'): return "TVN_"
    if term.startswith('Colorado deer associated narnavirus'): return 'CDaN_'
    if term.startswith('Culicoides Partiti-like'): return "CPlV_"
    if term.find('Bovine hepacivirus') > -1: return "BHC_"
    raise ValueError
    return "UNKNOWN_"

def parse_seq_header(header_line:str, gene_name:str) -> dict:
    gene = gene_name.replace('CDS','').strip()
    if gene.startswith('RNA-dependent'):
        gene = 'RdRP'
    elif gene.startswith('Hypothetical protein'):
        gene = ""
    elif gene.startswith('Putative capsid protein'):
        gene = ""
        
    fields = header_line.split('isolate')
    separator = gene_name.split(' ')[0]
    split_fields = fields[1].split(separator)
    virus_prefix = get_virus_prefix(fields[0])
    canon = virus_prefix + split_fields[0].strip().replace(" ", "_").replace('/','_')
    if canon.find(gene) == -1:
        canon = canon + '_' + gene.replace(' ', '_')
    return { 'canonical': canon }

def print_parsed_header(parsed_dict):
    for k,v in parsed_dict.items():
        print(f"{k}={v}", end=' ')

    print()


def convert_date_and_location(row:dict) -> dict:

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

def validate_genbank_fields(fields:dict) -> dict:
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

import os, glob
def read_all_Gffs() -> pd.DataFrame:
    return pd.concat([read_Aim2_Gffs(), read_Aim3_Gffs()], axis = 0)

def read_Aim2_Gffs() -> pd.DataFrame:
    gff = pd.DataFrame()
    for path in glob.glob(f"{os.path.dirname(sys.argv[0])}/**/Aim_2_SRA_Stuff/**/*_*.gff", recursive=True):
        this_gff = read_GFF(path)
        this_gff['source_file'] = path
        gff = pd.concat([gff, this_gff], axis = 0)

    return gff
def read_Aim3_Gffs() -> pd.DataFrame:
    gff = pd.DataFrame()
    for path in glob.glob(f"{os.path.dirname(sys.argv[0])}/**/Aim_3_SRA_Stuff/**/*_*.gff", recursive=True):
        this_gff = read_GFF(path)
        this_gff['source_file'] = path
        gff = pd.concat([gff, this_gff], axis = 0)

    return gff

def read_Aim2_Fastas() -> dict:
    all_seqs = {}
    for path in glob.glob(f"{os.path.dirname(sys.argv[0])}/**/Aim_2_SRA_Stuff/**/*_*.fasta", recursive=True):
        for i, record in enumerate(SeqIO.parse(path, "fasta")):
            seq_id = record.description.strip().replace(' ','_').replace('/', '_').replace('Bluetongue_virus_isolate_','')
            all_seqs[seq_id] = str(record.seq)
            # if seq_id in metadata:
            #     print('>' + seq_id.replace(' ', '_'), end=' ')
            #     for h in header:
            #         print(f"[{h}={metadata[seq_id][h]}]", end=" ")
            #     print()
            # else:
            #     print(f"Warning: {seq_id} not in data file", file=sys.stderr)
            #     continue

            #print(str(record.seq))

    return all_seqs

def read_Aim3_Fastas() -> dict:
    all_seqs = {}
    for path in glob.glob(f"{os.path.dirname(sys.argv[0])}/**/Aim_3_SRA_Stuff/**/*_*.fasta", recursive=True):
        for i, record in enumerate(SeqIO.parse(path, "fasta")):
            seq_id = record.description.strip()
            all_seqs[seq_id] = str(record.seq)
            # if seq_id in metadata:
            #     print('>' + seq_id.replace(' ', '_'), end=' ')
            #     for h in header:
            #         print(f"[{h}={metadata[seq_id][h]}]", end=" ")
            #     print()
            # else:
            #     print(f"Warning: {seq_id} not in data file", file=sys.stderr)
            #     continue

            print(str(record.seq))

    return all_seqs

def main():
    # for testing only
    
    for path in glob.glob(f"{os.path.dirname(sys.argv[0])}/**/Aim_*_SRA_Stuff/**/*_*.gff", recursive=True):
        gff = read_GFF(path)
        gff.to_csv(sys.stdout, sep="\t", index=False,header=False)


if __name__ == "__main__": main()
