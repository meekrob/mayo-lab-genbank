#!/usr/bin/env python3
import sys
import os
import glob
from urllib.parse import unquote # Geneious has escaped the spaces in URL encoding

with open("concatenated.gff", "w") as outfile:
    seen_headers = set()

    for filepath in glob.glob("GFF_Annotations/**/*.gff"):
        print(f"Processing {filepath}", file=sys.stderr)
        with open(filepath, "r") as infile:
            for line in infile:
                decoded = unquote(line)
                
                # Write headers only once
                if decoded.startswith("##"):
                    if decoded not in seen_headers:
                        outfile.write(decoded)
                        seen_headers.add(decoded)
                else:
                    # Always write data lines
                    fields = decoded.split("\t", 1)  # Split into max 2 parts
                    if len(fields) == 2:
                        seq_id = fields[0].replace(" ", "_")
                        outfile.write(seq_id + "\t" + fields[1])
                    else:
                        outfile.write(decoded)
