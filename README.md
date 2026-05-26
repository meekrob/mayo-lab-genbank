# Scripts for genbank submissions

The following code has custom scripts that convert sequence names to a Genbank-ready 
format in three different file types: GFF (exported from Geneious), FASTA (exported from Geneious), TSV/Excel (user's annotations).

The submitters are Tillie Dunham and Mollie Burton, whose sequences and annotation formats are similar but require individual attention.

```mermaid
graph LR;
    Py_in@{ shape: tri, label: "map.py\nin" }
    Py_out@{ shape: manual-file, label: "map.py\nout" }
    user_annot@{ shape: manual-input, label: "TSV/Excel" }
    Fasta@{ shape: docs }
    GFF@{ shape: docs }
    user_variation_geneious@{ shape: subproc, label: "Detect user label conventions in Geneious"}
    user_variation_xl@{ shape: subproc, label: "Detect user label conventions in Excel Spreadsheet"}
    canonical_seqid@{ shape: subproc, label: "Genbank-ready seq IDs"}

    CDS@{ shape: bow-rect, label: "CDS ranges" }
    sequences@{ shape: docs }
    metadata@{ shape: bow-rect, label: "metadata"}; 


    GFF-->map_tsv_to_fasta_header.py;
    user_annot-->map_tsv_to_fasta_header.py;
    Fasta-->map_tsv_to_fasta_header.py;

    subgraph map_tsv_to_fasta_header.py
        direction LR
        Py_in-->seq_parser;
        Py_in-->gff_parser;
        Py_in-->metadata_parser;
        seq_parser-->user_variation_geneious;
        seq_parser-->sequences;
        gff_parser-->user_variation_geneious;
        gff_parser-->CDS;
        canonical_seqid-->CDS;
        metadata_parser-->user_variation_xl;
        metadata_parser-->metadata;
        canonical_seqid-->metadata;
        canonical_seqid-->sequences;
        user_variation_xl-->canonical_seqid;
        user_variation_geneious-->canonical_seqid;
        Py_in== User-specific data module ==>Py_out;
        sequences-->validate_seq;
        sequences-->add_modifiers;
        CDS-->validate_seq;
        metadata-->add_modifiers;
        
    end
    validate_seq-->tbl;
    add_modifiers-->fsa;
    table2asn-->sqn;
    tbl-->table2asn;
    fsa-->table2asn;
    
    table2asn-->val;

```
