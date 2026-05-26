# Scripts for genbank submissions

The following code has custom scripts that convert sequence names to a Genbank-ready 
format in three different file types: GFF (exported from Geneious), FASTA (exported from Geneious), TSV/Excel (user's annotations).

The submitters are Tillie Dunham and Mollie Burton, whose sequences and annotation formats are similar but require individual attention.

```mermaid
graph LR;
    user_annot@{ shape: manual-input, label: "TSV/Excel" }
    Fasta@{ shape: docs }
    GFF@{ shape: docs }
    user_variation_geneious@{ shape: subproc, label: "Detect user label conventions in Geneious"}
    user_variation_xl@{ shape: subproc, label: "Detect user label conventions in Excel Spreadsheet"}
    canonical_seqid@{ shape: subproc, label: "canonicalize seqid"}
    add_modifiers@{ shape: subproc}
    validate_seq@{ shape: subproc}
    table2asn@{ shape: diamond}

    CDS@{ shape: bow-rect, label: "CDS ranges" }
    sequences@{ shape: docs }
    metadata@{ shape: bow-rect, label: "metadata"}; 

    subgraph "user data"
    GFF;
    user_annot;
    Fasta;
    end

    
    GFF-->gff_parser;
    user_annot-->metadata_parser;
    Fasta-->seq_parser;


    subgraph scripts
        direction LR
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

    subgraph "submission files"
        table2asn
        tbl
        fsa
        sqn
    end

```
