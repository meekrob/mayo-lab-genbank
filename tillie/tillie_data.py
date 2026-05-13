# tillie_data.py
import sys
import re
from datetime import datetime
from genbank_vals import geoloc_countries, state_names_to_abbrev, host_lookup
import pandas as pd
from urllib.parse import unquote # Geneious output is URL-quoted

def get_gene_name(attributes) -> str | None:
    match = re.search(r'Name=(\S+)', attributes)
    if match and match.group(1).strip():
        return match.group(1)
    return None

def get_seq_id_prot(orig_seq_id:str) -> str:
    # this is a special function to get just the gene name as specified by the seq ID in the GFF file. 
    # It needs to be called separately in order to contrast it with the gene_name parsed from the attribute column
    parsed = parse_seq_header(orig_seq_id)
    return parsed['prot']

def read_GFF(gff_file) -> pd.DataFrame:
    # gff
    GFF = pd.read_csv(gff_file, sep="\t", header = None, comment='#')
    GFF.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    # processing and filtering
    GFF = GFF[ (GFF['type'] == 'CDS') & 
              (~GFF['attributes'].str.contains(r'Name=$', regex=True))].copy() # filter out Geneious "Editing History..." feature rows
                                                                               # and empty Name= attribute rows.

    # extract gene name from attributes, but fail if not found
    GFF['gene_name'] = GFF['attributes'].apply(get_gene_name)
    missing_gene_name = GFF[GFF['gene_name'].isna()]
    if not missing_gene_name.empty:
        print("CDS rows with unparseable gene names:", file=sys.stderr)
        print(missing_gene_name[['seqid', 'attributes']], file=sys.stderr)
        raise ValueError(f"{len(missing_gene_name)} CDS rows failed gene name extraction")
    
    # parse seqid into our canonical one
    GFF['orig_id'] = GFF['seqid']
    GFF['seqid'] = GFF['seqid'].apply(unquote)
    GFF['seqid'] = GFF.apply(
        lambda row: canonical_seq_name(row['seqid'], row['gene_name']),
        axis = 1
    )

    # identify the misnamed entries and correct them
    # GFF['prot'] = GFF['seqid'].apply(get_seq_id_prot)
    # misnamed_rows = GFF['prot'] != GFF['gene_name']
    # if len(misnamed_rows) > 0:
    #     raise ValueError(f"There are still {len(misnamed_rows)}")


    # multi-annotated CDS rows (only 56, but still an issue), differ by length. Record span to take the longer one
    GFF['span'] = GFF['end'] - GFF['start']
    GFF = GFF.sort_values('span', ascending=False)

    # Log the discarded duplicates before dropping them
    discarded = GFF[GFF.duplicated(subset=['seqid', 'gene_name'], keep='first')]
    if not discarded.empty:
        print(f"Discarding {len(discarded)} redundant CDS annotations; keeping the longest entry.", file=sys.stderr)
        #print(discarded[['seqid', 'gene_name', 'start', 'end', 'span']], file=sys.stderr)

    GFF = GFF.drop_duplicates(subset=['seqid', 'gene_name'], keep='first')
    GFF = GFF.drop(columns='span')

    # check for duplicate records and crash if found
    dupes = GFF.groupby(['seqid', 'gene_name']).size()
    dupes = dupes[dupes > 1]
    if not dupes.empty:
        raise ValueError(f"Genuine duplicate (seqid, gene_name) pairs found:\n{dupes}")
    
    return GFF

def canonical_seq_name_from_parsed(parsed:dict):
    if 'date' in parsed:
        date_str = parsed['date'].replace('-','_').replace('.','_') + "_"
    else:
        date_str = ""

    if 'gene_name' in parsed:
        if 'prot' in parsed and parsed['prot'] != parsed['gene_name']:
            print(f"Info: correcting prot '{parsed['prot']}' -> '{parsed['gene_name']}' via Name attribute", file=sys.stderr)
    
        parsed['prot'] = parsed['gene_name'] # this dijection comes from GFFs where we had to get 
                                                                   # the gene name out of the attribute field
                                                                   # Caused by: overlapping CDSs

    canonical_name = f"Bluetongue_virus_{parsed['strain']}_{date_str}{parsed['prot']}"
    return canonical_name

def canonical_seq_name(header_line:str, gene_name:str | None = None):
    parsed = parse_seq_header(header_line)
    if gene_name is not None:
        parsed['prot'] = gene_name
    return canonical_seq_name_from_parsed(parsed)

def canonical_host(host_field):
    host = host_field.lower()
    if host in host_lookup:
        return host_lookup[host]
    
    return host

def id_seqid_type(fields, btv_col) -> dict:
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
            'strain': strain.replace('_','') # needed in sequence headers
        }
    
    # VP3     Clinical_16     BTV6    mule_deer       5.5yo   female  2021.10.19      2/7
    if btv_col == 2 and len(fields) == 7:
        prot, strain, btv, host, age, sex, date = fields
        return { 
            'prot': prot,
            'strain': strain,
            'btv': btv,
            'host': canonical_host(host),
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
            'host': canonical_host(fields[4]),
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
    
    if btv_col == 1 and len(fields) == 6:
        if(fields[2].startswith('A') or fields[2] == 'NA'): # this version is in the metadata file
            prot, btv, genotype, strain, host, date = fields
            return {
                'prot': prot,
                'btv': btv,
                'genotype': genotype,
                'strain': strain,
                'host': canonical_host(host),
                'date': date
            }
    
        else: # btv_col == 1 and len(fields) == 6: # this version is in the sequence header and GFF
            prot, btv, strain, host, location, date = fields
            return {
                'prot': prot,
                'btv': btv,
                'strain': strain,
                'location': location,
                'date': date
                
            }
    
    
    raise ValueError
    return { 'error': ';'.join(fields)}

def parse_seq_header(header_line:str) -> dict:
    fields = [field for field in re.split('[/|]', header_line.strip()) if field != '']
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
