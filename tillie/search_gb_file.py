#!/usr/bin/env python3
"""
search_gb_file.py
Parse a genbank file, extract selected IDs
"""
import sys
from Bio import SeqIO,BiopythonParserWarning
import warnings
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL) # omits broken pipe if piping to head and ctrl-c

progname = sys.argv[0]
USAGE = f"""{progname} - 
 Usage 1: Search genbank file by id:
    {progname} file.gb ID[,...]

 Usage 2: Print all ids from a genbank file
    {progname} file.gb 
"""

bluetongue_lineage = [
    "Viruses", 
    "Riboviria", 
    "Orthornavirae", 
    "Duplornaviridae", 
    "Sedoreovirinae", 
    "Orbivirus"
]


if len(sys.argv) < 2:
    print(USAGE)
    sys.exit(0)

gb_file = sys.argv[1]
IDs = []
if len(sys.argv) > 2:
    IDs = sys.argv[2:]

all_gb = {}
warnings.simplefilter('ignore', BiopythonParserWarning) # apparently now the LOCUS is wrong ... but it doesn't say why


for i, record in enumerate(SeqIO.parse(gb_file, "genbank")):
    # add some missing data
    if 'data_file_division' not in record.annotations: record.annotations["data_file_division"] = "VRL"
    if 'molecular_annotations' not in record.annotations: record.annotations['molecule_type'] = "RNA"
    if 'taxonomy' not in record.annotations: record.annotations["taxonomy"] = bluetongue_lineage
    all_gb[record.id] = record

    if IDs:
        if record.id in IDs:
            SeqIO.write(record, sys.stdout, "genbank")
    else:
            try:
                print(record.id)
            except BrokenPipeError:
                break