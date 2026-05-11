#!/usr/bin/env python3
import sys
import re

def id_seqid_type(fields, btv_col):
    if btv_col == -1:
        prot, strain = fields
        return {
            'prot': prot,
            'strain': strain
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

def print_parsed_header(parsed_dict):
    for k,v in parsed_dict.items():
        print(f"{k}=\"{v}\"", end=' ')

    print()

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

host_terms_to_sp = {
    'mule_deer': 'mule deer',
    'white_tailed_deer': 'white-tailed deer',
    'cow': 'cow',
    'bovine': 'cow',
    'sheep': 'sheep',
    'llama': 'llama',
    'bighorn_sheep': 'bighorn sheep',
    'reindeer': 'reindeer'
}

with open('/Users/david/work/mayo_lab_sequence_submission/tillie/new_files/TD_n1257_GenBank.fasta') as infh:
    for line in infh:
        if line.startswith('>'):
            parsed = parse_seq_header(line.lstrip('>'))
            if 'host' in parsed:
                if parsed['host'] in host_terms_to_sp: 
                    parsed['host'] = host_terms_to_sp[ parsed['host'] ]
                else:
                    del parsed['host']

            print_parsed_header(parsed)
            
            
