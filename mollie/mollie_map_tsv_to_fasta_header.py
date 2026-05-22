#!/usr/bin/env python3
import sys,os
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
from typing import TextIO
import mollie_data as mollie

### hard-code paths- this is not a generalizable script
basedir= os.path.dirname(sys.argv[0]) # where are we?
## input metadata file paths
aim2_tsv_file = f"{basedir}/SRA_materials/Aim_2_SRA_Stuff/Aim_2_SRA_metadata.xlsx"
aim3_tsv_file = f"{basedir}/SRA_materials/Aim_3_SRA_Stuff/Aim_3_SRA_metadata.xlsx"
## output file paths
### Aim 2
aim2_tbl_file = f"{basedir}/Aim2_submission_files/aim2_feature.tbl"
aim2_fsa_file = f"{basedir}/Aim2_submission_files/aim2_feature.fsa"
### Aim 3
aim3_tbl_file = f"{basedir}/Aim3_submission_files/aim3_feature.tbl"
aim3_fsa_file = f"{basedir}/Aim3_submission_files/aim3_feature.fsa"

def dfprint(df, out=sys.stderr, sep ="\t"):
    df.to_csv(out,sep)

def write_tbl_entry(genbank_id:str, cds_entries:list, tbl_out: TextIO) -> None:
    print(f">Feature {genbank_id}", file=tbl_out)
    for cds_entry in cds_entries:
        partial3prime: bool = cds_entry['partial3prime']
        partial5prime: bool = cds_entry['partial5prime']
        start: int = cds_entry['start']
        end: int = cds_entry['end']
        codon_start: int = cds_entry['codon_start']
        gene_name:str = cds_entry['gene_name']
        s = f"<{start}" if partial5prime else str(start)
        e = f">{end}"   if partial3prime else str(end)
    
        print(s, e, "CDS", sep="\t", file=tbl_out)
        print("\t" * 3 + "product", gene_name, sep="\t", file=tbl_out)
        print("\t" * 3 + "codon_start", codon_start, sep="\t", file=tbl_out)
        #print("\t" * 3 + "transl_table", 1, sep="\t", file=tbl_out)

def examine_sequence(seq:str, start:int, end:int, phase:int) -> tuple:
    if end > len(seq):
        raise ValueError
    cds = seq[(start-1)+phase:end]  # convert to 0-based
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    has_start = cds[:3].upper() == 'ATG'
    has_stop  = cds[-3:].upper() in stop_codons
    is_5prime_partial = not has_start
    is_3prime_partial = not has_stop
    
    codon_start = phase + 1
    return is_5prime_partial, is_3prime_partial, codon_start

def make_tbl_entry(cds_row, seq):
    start = int(cds_row['start'])
    end = int(cds_row['end'])
    gene_name = cds_row['gene_name']
    phase = frame_lookup[cds_row['phase']]
    # match annotation to sequence
    is_5prime_partial, is_3prime_partial, codon_start = examine_sequence(seq, start, end, phase)

    return { 
        'start': start,
        'end': end,
        'gene_name': gene_name,
        'partial5prime': is_5prime_partial, 
        'partial3prime': is_3prime_partial, 
        'codon_start': codon_start}

def write_fsa_entry(fasta_id, metarow, fsa_fh):

        print(f">{fasta_id}", file=fsa_fh, end=' ')

        annotations = []
        for k,v in metarow.items():
            if k in export_fields_keys:
                if v.item() == 'NA': continue

                modifier = f"[{k}={v.item()}]"
                annotations.append(modifier)

        print(" ".join(annotations), file=fsa_fh)
        print(seq, file=fsa_fh)

### global fields and data
frame_lookup = { '.':0, 
                '0': 0, 
                '1': 1, 
                '2': 2,
                 0: 0, 
                 1: 1, 
                 2: 2
                }

platform_mapping = {
    "Oxford Nanopore Technologies": "nanopore",
    "Illumina Sequencing": "illumina" # not in Mollie's data
}

export_fields_keys = ["genotype", "host", "strain", 
                      "collected_by", "geo_loc_name", "collection_date", 
                      "note", "isolation_source", "organism"] 

metadata = defaultdict(lambda: {}) # map by seq_id: { header_field: value }

### common munging operations needed on boths metadata sheets
### most of them are vectorized pandas.DataFrame operations
### put them together temporarily so the same code can be used while keeping modifications in-place
dfs = {
    'aim2' : pd.read_excel(aim2_tsv_file).dropna(how='all'),
    'aim3' : pd.read_excel(aim3_tsv_file).dropna(how='all')
}

for name, df in dfs.items():
    df['gene_name'] = df['sample_ID'].apply(mollie.get_gene_name_from_spreadsheet_sampleID)    
    df['seq_fa'] = df['consensus_sequence.fasta'].apply(lambda x: x.replace(' ', '_').replace('/','_').replace('.fasta',''))
    df['collection_date'] = df['collection_date'].dt.strftime("%d-%b-%Y")
    df.rename(columns={'source': 'isolation_source'}, inplace=True)

    # sequencing_platform isn't a modifier, but you can put it in a "note"
    df['note'] = df['sequencing_platform'].apply(lambda x: f"sequenced with {platform_mapping[x]}")
    df.drop(['sequencing_platform'], axis=1, inplace=True)
    
    df['seqid'] = df.apply(
    lambda row: mollie.canonical_seq_name(row['sample_ID'], row['gene_name']),
    axis=1
)
# separate back out
aim2 = dfs['aim2']
aim3 = dfs['aim3']

aim2['organism'] = "Bluetongue virus"

### Aim 2
if False:
    aim2_fa = mollie.read_Aim2_Fastas()
    gff2 = mollie.read_Aim2_Gffs()

    with open(aim2_tbl_file,"w") as aim2_tbl, open(aim2_fsa_file, "w") as aim2_fsa:

        for source_file, group in gff2.groupby('source_file'):
            metarow = aim2[aim2['seqid'] == group.iloc[0]['seqid'] ]
            fasta_id = group.iloc[0]['seqid']
            
            seq = aim2_fa[fasta_id]
            tbl_entries = []
            for idx, cds_row in group.iterrows():
                tbl_entries.append(make_tbl_entry(cds_row, seq))

            write_tbl_entry(fasta_id, tbl_entries, aim2_tbl)
            write_fsa_entry(fasta_id, metarow, aim2_fsa)

    print(f"Wrote:\n\tfsa={aim2_fsa_file},\n\ttbl={aim2_tbl_file}", file=sys.stderr)

### Aim 3
aim3_fa = mollie.read_Aim3_Fastas()
gff3 = mollie.read_Aim3_Gffs()

with open(aim3_tbl_file,"w") as aim3_tbl, open(aim3_fsa_file, "w") as aim3_fsa:

    for source_file, group in gff3.groupby('source_file'):
        metarow = aim3[aim3['seqid'] == group.iloc[0]['seqid'] ]
        fasta_id = group.iloc[0]['seqid']
        
        seq = aim3_fa[fasta_id]
        tbl_entries = []
        for idx, cds_row in group.iterrows():
            tbl_entries.append(make_tbl_entry(cds_row, seq))
        if fasta_id == 'TVN_CO_USA_2021_pol':
            # manual case of ribosome slipping that doesn't register as a polyprotein in the usual workflow
            TAB="\t"
            print(f""">Feature TVN_CO_USA_2021_pol
16	1938	gene
1938	4526	
			gene	pol
16	1938	CDS
1938	4526	
			gene	pol
			product	polyprotein
			exception	ribosomal slippage
			codon_start	1
16	1938	gene
			gene	cap
16	1938	CDS
			gene	cap
			product	capsid protein
			codon_start	1
16	4526	gene
			gene	pol
<1938	4526	CDS
			gene	pol
			product	RNA-dependent RNA polymerase
			codon_start	1
""", file=aim3_tbl)
        else:
            write_tbl_entry(fasta_id, tbl_entries, aim3_tbl)
        write_fsa_entry(fasta_id, metarow, aim3_fsa)

print(f"Wrote:\n\tfsa={aim3_fsa_file},\n\ttbl={aim3_tbl_file}", file=sys.stderr)
