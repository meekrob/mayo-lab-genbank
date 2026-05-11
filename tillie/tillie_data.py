# tillie_data.py
import sys
import re
from datetime import datetime
from genbank_vals import geoloc_countries, state_names_to_abbrev, host_lookup


def canonical_seq_name(header_line):
    parsed = parse_seq_header(header_line)
    canonical_name = f"Bluetongue_virus_{parsed['strain']}_{parsed['prot']}"
    return canonical_name

def canonical_host(host_field):
    host = host_field.lower()
    if host in host_lookup:
        return host_lookup[host]
    
    return host

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
