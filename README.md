# Scripts for genbank submissions

This project has custom scripts that convert sequence names to a Genbank-ready 
format in three different file types: GFF (exported from Geneious), FASTA (exported from Geneious), TSV/Excel (user's annotations),
and prepares files as input for *table2asn*.

The submitters are Tillie Dunham and Mollie Burton, whose sequences and annotation formats are similar but require individual attention.

```mermaid
graph LR;
    classDef validation fill:#f9f,stroke:#333,stroke-width:4px;
    class validate_seq,canonical_seqid,validate_host,validate_geo-loc,validate_date_format validation;
    user_annot@{ shape: manual-input, label: "TSV/Excel" }
    Fasta@{ shape: docs }
    GFF@{ shape: docs }
    user_variation_geneious@{ shape: subproc, label: "Detect user's conventions (Geneious)"}
    user_variation_xl@{ shape: subproc, label: "Detect user's conventions (Excel Spreadsheet)"}
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

        validate_host([validate host])
        validate_date_format([validate date])
        validate_geo-loc([validate geo-loc])

        metadata<-->validate_date_format;

        metadata<-->validate_geo-loc;

        metadata<-->validate_host;

        A@{ shape: sm-circ, label: "Small start" }
        validate_host-->A
        validate_geo-loc-->A
        validate_date_format-->A

        canonical_seqid-->metadata;
        canonical_seqid-->sequences;
        user_variation_xl-->canonical_seqid;
        user_variation_geneious-->canonical_seqid;
        sequences-->validate_seq;
        sequences-->add_modifiers;
        CDS-->validate_seq;
        A-->add_modifiers;
        
    end
    A-->tsv;
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

    table2asn-->gbn;
    table2asn-->val;

    subgraph "extra output"
        direction LR
        tsv;
        gbn;
        val;
    end



```
